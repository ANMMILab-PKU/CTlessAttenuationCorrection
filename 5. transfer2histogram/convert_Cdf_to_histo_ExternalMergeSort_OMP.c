#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h> // 引入 OpenMP 头文件

// ================= 配置区域 =================
// 每个线程处理的块大小 (1000万事件 ~ 120MB)
#define CHUNK_SIZE 10000000 
#define TEMP_FILE_PREFIX "temp_"

// ================= 数据结构 (保持不变) =================
#pragma pack(push, 1)
typedef struct {
    uint32_t time;
    uint32_t id1;
    uint32_t id2;
} RawEvent;

typedef struct {
    uint32_t counts;
    uint32_t id1;
    uint32_t id2;
} HistoEvent;
#pragma pack(pop)

// ================= 辅助函数 (保持不变) =================
int compare_events(const void *a, const void *b) {
    HistoEvent *ea = (HistoEvent *)a;
    HistoEvent *eb = (HistoEvent *)b;
    if (ea->id1 != eb->id1) return (ea->id1 < eb->id1) ? -1 : 1;
    if (ea->id2 != eb->id2) return (ea->id2 < eb->id2) ? -1 : 1;
    return 0;
}

void normalize_ids(uint32_t *id1, uint32_t *id2) {
    if (*id1 > *id2) {
        uint32_t temp = *id1;
        *id1 = *id2;
        *id2 = temp;
    }
}

// ================= 并行阶段 1: 多线程分块处理 =================

int process_chunks_parallel(const char *input_path, long long *total_events_read) {
    FILE *fin = fopen(input_path, "rb");
    if (!fin) { perror("Error opening input file"); exit(1); }

    // 获取系统最大线程数，例如 8核 CPU 通常设为 8 或 16
    int max_threads = omp_get_max_threads();
    printf("--- Phase 1: Parallel Processing with %d Threads ---\n", max_threads);

    // 分配一个巨大的缓冲区，一次性读入足够所有线程处理的数据
    // 总内存占用 = 线程数 * CHUNK_SIZE * sizeof(RawEvent)
    size_t batch_buffer_size = (size_t)max_threads * CHUNK_SIZE;
    RawEvent *raw_batch_buffer = (RawEvent *)malloc(sizeof(RawEvent) * batch_buffer_size);
    
    // 输出缓冲区，每个线程只需要访问自己的部分
    HistoEvent *sort_batch_buffer = (HistoEvent *)malloc(sizeof(HistoEvent) * batch_buffer_size);

    if (!raw_batch_buffer || !sort_batch_buffer) {
        fprintf(stderr, "Memory allocation failed! Reduce CHUNK_SIZE.\n");
        exit(1);
    }

    int batch_idx = 0;
    int global_chunk_count = 0; // 记录总共生成了多少个临时文件
    *total_events_read = 0;

    while (1) {
        // 主线程负责读取一大批数据 (IO 操作无法并行，单线程读最快)
        size_t items_read = fread(raw_batch_buffer, sizeof(RawEvent), batch_buffer_size, fin);
        if (items_read == 0) break;

        *total_events_read += items_read;
        
        // 计算这一批数据能分多少个完整的 chunk，以及是否有多余的零头
        // 为了简化并行逻辑，我们将读取到的 items_read 均匀分配给各个线程
        // 注意：这里我们按线程切分任务，每个线程处理一部分数组
        
        int chunks_in_batch = max_threads; 
        
        // 开启并行区域
        #pragma omp parallel for
        for (int t = 0; t < chunks_in_batch; t++) {
            // 计算当前线程负责的数据范围
            size_t start_idx = t * (items_read / chunks_in_batch);
            size_t end_idx = (t + 1) * (items_read / chunks_in_batch);
            
            // 最后一个线程负责处理剩下的余数 (如果有)
            if (t == chunks_in_batch - 1) {
                end_idx = items_read;
            }

            size_t local_count = end_idx - start_idx;
            if (local_count > 0) {
                // 指针定位到当前线程的缓冲区位置
                RawEvent *local_raw = raw_batch_buffer + start_idx;
                HistoEvent *local_sort = sort_batch_buffer + start_idx;

                // 1. 标准化 (Thread Local)
                for (size_t i = 0; i < local_count; i++) {
                    uint32_t u = local_raw[i].id1;
                    uint32_t v = local_raw[i].id2;
                    normalize_ids(&u, &v);
                    local_sort[i].counts = 1;
                    local_sort[i].id1 = u;
                    local_sort[i].id2 = v;
                }

                // 2. 排序 (Thread Local - 并行加速的核心)
                // 各个线程同时对自己的一小段内存进行快排
                qsort(local_sort, local_count, sizeof(HistoEvent), compare_events);

                // 3. 压缩并写入临时文件
                // 文件名需要包含 batch 和 thread id 以避免冲突
                char temp_filename[256];
                int unique_id = batch_idx * max_threads + t;
                sprintf(temp_filename, "%s%d.tmp", TEMP_FILE_PREFIX, unique_id);
                
                FILE *ftemp = fopen(temp_filename, "wb");
                if (ftemp) {
                    uint32_t cur_id1 = local_sort[0].id1;
                    uint32_t cur_id2 = local_sort[0].id2;
                    uint32_t cur_cnt = local_sort[0].counts;

                    for (size_t i = 1; i < local_count; i++) {
                        if (local_sort[i].id1 == cur_id1 && local_sort[i].id2 == cur_id2) {
                            cur_cnt += local_sort[i].counts;
                        } else {
                            HistoEvent out = {cur_cnt, cur_id1, cur_id2};
                            fwrite(&out, sizeof(HistoEvent), 1, ftemp);
                            cur_id1 = local_sort[i].id1;
                            cur_id2 = local_sort[i].id2;
                            cur_cnt = local_sort[i].counts;
                        }
                    }
                    HistoEvent out = {cur_cnt, cur_id1, cur_id2};
                    fwrite(&out, sizeof(HistoEvent), 1, ftemp);
                    fclose(ftemp);
                } else {
                    fprintf(stderr, "Error creating temp file %s\n", temp_filename);
                }
            }
        } // 隐式同步：主线程会等待所有子线程完成

        // 更新总块数
        global_chunk_count += chunks_in_batch;
        printf("Batch %d processed. (Total events so far: %lld)\n", batch_idx, *total_events_read);
        batch_idx++;
    }

    free(raw_batch_buffer);
    free(sort_batch_buffer);
    fclose(fin);

    // 返回生成的最大可能的文件索引数（注意：如果最后一次读取数据很少，可能部分线程没生成文件，需要在merge时判断文件是否存在）
    return batch_idx * max_threads;
}

// ================= 阶段 2: 归并 (IO密集型，保持单线程) =================
// 注意：需要微调 merge 函数以处理可能不存在的临时文件
void merge_chunks(int max_potential_chunks, const char *output_path) {
    printf("--- Phase 2: Merging (Serial) ---\n");

    // 动态分配指针数组
    FILE **temp_files = (FILE **)calloc(max_potential_chunks, sizeof(FILE*));
    HistoEvent *current_heads = (HistoEvent *)malloc(sizeof(HistoEvent) * max_potential_chunks);
    int *has_data = (int *)calloc(max_potential_chunks, sizeof(int));
    
    int active_files_count = 0;

    // 打开存在的临时文件
    for (int i = 0; i < max_potential_chunks; i++) {
        char temp_filename[256];
        sprintf(temp_filename, "%s%d.tmp", TEMP_FILE_PREFIX, i);
        temp_files[i] = fopen(temp_filename, "rb");
        
        if (temp_files[i]) {
            size_t r = fread(&current_heads[i], sizeof(HistoEvent), 1, temp_files[i]);
            if (r > 0) {
                has_data[i] = 1;
                active_files_count++;
            } else {
                fclose(temp_files[i]);
                temp_files[i] = NULL;
            }
        }
    }
    printf("Found %d valid chunks to merge.\n", active_files_count);

    FILE *fout = fopen(output_path, "wb");
    if (!fout) { perror("Error creating output file"); exit(1); }

    // 归并循环 (使用简单的线性查找最小值，如果 active_files_count > 100，建议改用最小堆 Min-Heap)
    while (active_files_count > 0) {
        int min_idx = -1;
        uint32_t min_id1 = 0xFFFFFFFF;
        uint32_t min_id2 = 0xFFFFFFFF;

        // 寻找最小值
        for (int i = 0; i < max_potential_chunks; i++) {
            if (has_data[i]) {
                int is_smaller = 0;
                if (min_idx == -1) is_smaller = 1;
                else if (current_heads[i].id1 < min_id1) is_smaller = 1;
                else if (current_heads[i].id1 == min_id1 && current_heads[i].id2 < min_id2) is_smaller = 1;

                if (is_smaller) {
                    min_idx = i;
                    min_id1 = current_heads[i].id1;
                    min_id2 = current_heads[i].id2;
                }
            }
        }

        if (min_idx == -1) break;

        uint32_t total_counts = 0;
        
        // 累加所有等于最小值的 ID 对
        for (int i = 0; i < max_potential_chunks; i++) {
            if (has_data[i] && current_heads[i].id1 == min_id1 && current_heads[i].id2 == min_id2) {
                total_counts += current_heads[i].counts;
                size_t r = fread(&current_heads[i], sizeof(HistoEvent), 1, temp_files[i]);
                if (r == 0) {
                    has_data[i] = 0;
                    active_files_count--;
                    fclose(temp_files[i]);
                    temp_files[i] = NULL;
                    // 删除临时文件
                    char temp_filename[256];
                    sprintf(temp_filename, "%s%d.tmp", TEMP_FILE_PREFIX, i);
                    remove(temp_filename);
                }
            }
        }
        
        HistoEvent out = {total_counts, min_id1, min_id2};
        fwrite(&out, sizeof(HistoEvent), 1, fout);
    }

    free(current_heads);
    free(has_data);
    free(temp_files);
    fclose(fout);
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Usage: %s <input.Cdf> <output.histo>\n", argv[0]);
        return 1;
    }
    long long total_events = 0;
    
    // 1. 并行处理生成临时文件
    int max_chunks = process_chunks_parallel(argv[1], &total_events);
    
    // 2. 归并
    merge_chunks(max_chunks, argv[2]);

    printf("Done! Processed %lld events in parallel.\n", total_events);
    return 0;
}
