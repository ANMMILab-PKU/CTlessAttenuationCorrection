/*
 * transfer_sortedDat_into_Cdf_timeDiff_transmission.c
 *
 * 功能：筛选 511 emission 真符合事件
 * 输入：.dat 和 .consist 文件
 * 输出：.Cdf (1, c1, c2) 和 .timeDiff (dt)
 * 策略：分块同步读取，流式处理
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

// 定义批处理大小（例如每次处理 65536 个事件），根据内存大小可调
#define BATCH_SIZE 65536

// 输入 .dat 文件的数据结构 (24 bytes)
typedef struct {
    uint32_t crystalID1;
    uint32_t crystalID2;
    float timeDiff;
    float energy1;
    float energy2;
    uint32_t type;
} DatEvent;

// 输入 .consist 文件的数据结构 (24 bytes)
typedef struct {
    uint32_t eid_consist;
    uint32_t sid1;
    uint32_t sid2;
    uint32_t cp1;
    uint32_t cp2;
    uint32_t bs_type;
} ConsistEvent;

// 输出 .Cdf 文件的数据结构 (12 bytes)
typedef struct {
    uint32_t flag; // 固定为 1
    uint32_t c1;
    uint32_t c2;
} CdfOutput;

// 辅助函数：替换文件扩展名
void replace_extension(const char *original, const char *new_ext, char *output) {
    strcpy(output, original);
    char *dot = strrchr(output, '.');
    if (dot) {
        *dot = '\0'; // 去掉旧后缀
    }
    strcat(output, new_ext); // 加上新后缀
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <input_filename.dat>\n", argv[0]);
        fprintf(stderr, "Note: The program expects a corresponding .consist file to exist.\n");
        return EXIT_FAILURE;
    }

    char *input_dat_path = argv[1];
    char input_consist_path[512];
    char output_cdf_path[512];
    char output_td_path[512];

    // 1. 生成文件名
    replace_extension(input_dat_path, ".consist", input_consist_path);
    replace_extension(input_dat_path, ".Cdf", output_cdf_path);
    replace_extension(input_dat_path, ".timeDiff", output_td_path);

    printf("Input DAT:     %s\n", input_dat_path);
    printf("Input CONSIST: %s\n", input_consist_path);
    printf("Output CDF:    %s\n", output_cdf_path);
    printf("Output TD:     %s\n", output_td_path);

    // 2. 打开文件
    FILE *f_dat = fopen(input_dat_path, "rb");
    FILE *f_consist = fopen(input_consist_path, "rb");
    FILE *f_cdf = fopen(output_cdf_path, "wb");
    FILE *f_td = fopen(output_td_path, "wb");

    if (!f_dat || !f_consist || !f_cdf || !f_td) {
        perror("Error opening files");
        if (f_dat) fclose(f_dat);
        if (f_consist) fclose(f_consist);
        if (f_cdf) fclose(f_cdf);
        if (f_td) fclose(f_td);
        return EXIT_FAILURE;
    }

    // 3. 分配内存缓冲区
    DatEvent *buf_dat = (DatEvent *)malloc(BATCH_SIZE * sizeof(DatEvent));
    ConsistEvent *buf_consist = (ConsistEvent *)malloc(BATCH_SIZE * sizeof(ConsistEvent));
    
    // 输出缓冲（最大可能大小等于输入大小）
    CdfOutput *buf_out_cdf = (CdfOutput *)malloc(BATCH_SIZE * sizeof(CdfOutput));
    float *buf_out_td = (float *)malloc(BATCH_SIZE * sizeof(float));

    if (!buf_dat || !buf_consist || !buf_out_cdf || !buf_out_td) {
        perror("Memory allocation failed");
        return EXIT_FAILURE;
    }

    size_t total_processed = 0;
    size_t total_selected = 0;
    size_t n_read_dat, n_read_consist;

    printf("Processing...\n");

    // 4. 循环处理
    while (1) {
        // 同步读取两个文件
        n_read_dat = fread(buf_dat, sizeof(DatEvent), BATCH_SIZE, f_dat);
        n_read_consist = fread(buf_consist, sizeof(ConsistEvent), BATCH_SIZE, f_consist);

        // 校验文件对齐
        if (n_read_dat != n_read_consist) {
            fprintf(stderr, "\nError: File length mismatch! .dat has %zu items, .consist has %zu items in current batch.\n", 
                    n_read_dat, n_read_consist);
            break; // 或者选择处理较小的那个数量
        }

        if (n_read_dat == 0) break; // 处理完毕

        size_t batch_selected = 0;

        // 遍历当前 Batch 进行筛选
        for (size_t i = 0; i < n_read_dat; i++) {
            // 筛选条件：eid=1, sid1=0, sid2=0, cp1=0, cp2=0
            if(! (buf_consist[i].eid_consist == 1 &&
                buf_consist[i].sid1 == 0 &&
                buf_consist[i].sid2 == 0 &&
                buf_consist[i].cp1 == 0 &&
                buf_consist[i].cp2 == 0)) {
                
                // 填充 CDF 缓冲 (1, c1, c2)
                buf_out_cdf[batch_selected].flag = 1;
                buf_out_cdf[batch_selected].c1 = buf_dat[i].crystalID1;
                buf_out_cdf[batch_selected].c2 = buf_dat[i].crystalID2;

                // 填充 TimeDiff 缓冲
                buf_out_td[batch_selected] = buf_dat[i].timeDiff * 1.0e12f * (-1);

                batch_selected++;
            }
        }

        // 批量写入输出文件
        if (batch_selected > 0) {
            fwrite(buf_out_cdf, sizeof(CdfOutput), batch_selected, f_cdf);
            fwrite(buf_out_td, sizeof(float), batch_selected, f_td);
        }

        total_processed += n_read_dat;
        total_selected += batch_selected;
        
        // 简单的进度条效果
        // printf("\rProcessed: %zu, Selected: %zu", total_processed, total_selected);
    }

    printf("\nDone.\n");
    printf("Total Events Processed: %zu\n", total_processed);
    printf("Total Events Selected:  %zu (%.2f%%)\n", total_selected, 
           total_processed > 0 ? (double)total_selected / total_processed * 100.0 : 0.0);

    // 5. 清理资源
    free(buf_dat);
    free(buf_consist);
    free(buf_out_cdf);
    free(buf_out_td);
    fclose(f_dat);
    fclose(f_consist);
    fclose(f_cdf);
    fclose(f_td);

    return EXIT_SUCCESS;
}
