#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>

// ================= 配置区域 =================
#define EVENT_BUFFER_SIZE 4096  // 缓冲读取大小
#define MAX_FILES 10            // 最大支持文件数

// ================= 数据结构 =================
#pragma pack(push, 1)
typedef struct {
    uint32_t counts;
    uint32_t id1;
    uint32_t id2;
} HistoEvent;
#pragma pack(pop)

// 文件上下文
typedef struct {
    FILE *fin;
    FILE *fout;
    char out_filename[256];
    
    HistoEvent *buffer;
    size_t buf_count;
    size_t buf_pos;
    
    HistoEvent current_evt;
    int finished;
} FileContext;

// ================= 辅助函数 =================

// 比较函数：ID1 为主序，ID2 为次序
int compare_ids(const HistoEvent *a, const HistoEvent *b) {
    if (a->id1 < b->id1) return -1;
    if (a->id1 > b->id1) return 1;
    
    if (a->id2 < b->id2) return -1;
    if (a->id2 > b->id2) return 1;
    
    return 0;
}

// 预读函数
void fetch_next(FileContext *ctx) {
    if (ctx->finished) return;

    if (ctx->buf_pos >= ctx->buf_count) {
        size_t r = fread(ctx->buffer, sizeof(HistoEvent), EVENT_BUFFER_SIZE, ctx->fin);
        if (r == 0) {
            ctx->finished = 1;
            ctx->current_evt.id1 = UINT32_MAX;
            ctx->current_evt.id2 = UINT32_MAX;
            ctx->current_evt.counts = 0;
            return;
        }
        ctx->buf_count = r;
        ctx->buf_pos = 0;
    }
    ctx->current_evt = ctx->buffer[ctx->buf_pos++];
}

// 生成文件名
void create_output_filename(const char *input, char *output) {
    const char *ext = strrchr(input, '.');
    if (!ext) {
        sprintf(output, "%s_unified.histo", input);
    } else {
        size_t len = ext - input;
        strncpy(output, input, len);
        output[len] = '\0';
        strcat(output, "_unified.histo");
    }
}

// ================= 核心逻辑 =================

void unify_files(int num_files, char **filenames) {
    FileContext ctxs[MAX_FILES];
    
    printf("--- Processing %d files (Target: ID1 <= ID2) ---\n", num_files);

    // 1. 初始化
    for (int i = 0; i < num_files; i++) {
        ctxs[i].fin = fopen(filenames[i], "rb");
        if (!ctxs[i].fin) {
            fprintf(stderr, "Error opening input: %s\n", filenames[i]);
            exit(1);
        }

        create_output_filename(filenames[i], ctxs[i].out_filename);
        ctxs[i].fout = fopen(ctxs[i].out_filename, "wb");
        if (!ctxs[i].fout) {
            fprintf(stderr, "Error creating output: %s\n", ctxs[i].out_filename);
            exit(1);
        }

        ctxs[i].buffer = (HistoEvent *)malloc(sizeof(HistoEvent) * EVENT_BUFFER_SIZE);
        ctxs[i].buf_pos = 0;
        ctxs[i].buf_count = 0;
        ctxs[i].finished = 0;

        fetch_next(&ctxs[i]);
    }

    long long processed_pairs = 0;

    // 2. 循环归并
    while (1) {
        // 检查结束条件
        int all_finished = 1;
        for (int i = 0; i < num_files; i++) {
            if (!ctxs[i].finished) {
                all_finished = 0;
                break;
            }
        }
        if (all_finished) break;

        // A. 找当前最小值
        HistoEvent min_evt;
        min_evt.id1 = UINT32_MAX;
        min_evt.id2 = UINT32_MAX;

        for (int i = 0; i < num_files; i++) {
            if (!ctxs[i].finished) {
                if (compare_ids(&ctxs[i].current_evt, &min_evt) < 0) {
                    min_evt.id1 = ctxs[i].current_evt.id1;
                    min_evt.id2 = ctxs[i].current_evt.id2;
                }
            }
        }

        // B. 写入所有文件
        for (int i = 0; i < num_files; i++) {
            HistoEvent out_evt;
            
            // [修改处]：直接使用 min_evt 的 ID 顺序
            // 因为输入文件已经是 ID1 <= ID2，所以 min_evt 也是 ID1 <= ID2
            out_evt.id1 = min_evt.id1;
            out_evt.id2 = min_evt.id2;

            // 如果当前文件匹配到了这个 ID 组合
            if (!ctxs[i].finished && 
                ctxs[i].current_evt.id1 == min_evt.id1 && 
                ctxs[i].current_evt.id2 == min_evt.id2) {
                
                out_evt.counts = ctxs[i].current_evt.counts; // 写入真实计数
                fetch_next(&ctxs[i]); // 指针后移
            } else {
                out_evt.counts = 0; // 补零
            }

            fwrite(&out_evt, sizeof(HistoEvent), 1, ctxs[i].fout);
        }

        processed_pairs++;
        if (processed_pairs % 2000000 == 0) {
            printf("\rUnifying pairs: %lld ...", processed_pairs);
            fflush(stdout);
        }
    }

    printf("\nDone! Total unique pairs: %lld\n", processed_pairs);

    // 3. 清理
    for (int i = 0; i < num_files; i++) {
        fclose(ctxs[i].fin);
        fclose(ctxs[i].fout);
        free(ctxs[i].buffer);
        printf("Output: %s\n", ctxs[i].out_filename);
    }
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("Usage: %s file1.histo file2.histo [file3.histo ...]\n", argv[0]);
        return 1;
    }
    unify_files(argc - 1, &argv[1]);
    return 0;
}
