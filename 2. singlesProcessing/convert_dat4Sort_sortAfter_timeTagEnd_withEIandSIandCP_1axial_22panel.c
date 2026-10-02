#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

#define CHUNK_EVENTS 1000000
#define NB_CRYSTALS_PER_RING_L 1056
#define MAX_EN (1 << 9)
#define MAX_TM (1 << 26) // 更新为 26 位上限

// Helper: get file size in bytes
long get_file_size(const char *filename)
{
  struct stat st;
  if (stat(filename, &st) == 0)
    return st.st_size;
  return -1;
}

// Structure to hold one event for sorting
typedef struct
{
  uint64_t tmBits;
  uint32_t crID;
  float en;
  uint32_t evID;
  uint32_t soID;
  uint32_t cpID;
  size_t original_index;
} Event;

// Comparator for qsort
int compare_events(const void *a, const void *b)
{
  uint64_t t1 = ((Event *)a)->tmBits;
  uint64_t t2 = ((Event *)b)->tmBits;
  if (t1 < t2)
    return -1;
  if (t1 > t2)
    return 1;
  return 0;
}

void encode_singles(
    const char *crystalID_file,
    const char *time_file,
    const char *energy_file,
    const char *eventID_file,
    const char *sourceID_file,
    const char *comptonP_file,
    const char *output_dat,
    const char *refere_dat)
{
  const double ENERGY_OFFSET = 0.05;       // MeV
  const double ENERGY_RESOLUTION = 0.0025; // MeV/bin
  const double CFD_RESOLUTION = 7e-12;     // s

  // 更新为 26 位的 Time Tag 间隔阈值
  const uint64_t TIME_TAG_INTERVAL_CFD_RESOLUTION = (1ULL << 26);

  FILE *fid_cr = fopen(crystalID_file, "rb");
  FILE *fid_tm = fopen(time_file, "rb");
  FILE *fid_en = fopen(energy_file, "rb");
  FILE *fid_ev = fopen(eventID_file, "rb");
  FILE *fid_so = fopen(sourceID_file, "rb");
  FILE *fid_cp = fopen(comptonP_file, "rb");

  if (!fid_cr || !fid_tm || !fid_en || !fid_ev || !fid_so || !fid_cp)
  {
    fprintf(stderr, "Error opening input files.\n");
    exit(EXIT_FAILURE);
  }

  FILE *fid_out = fopen(output_dat, "wb");
  FILE *fid_ref = fopen(refere_dat, "wb");
  if (!fid_out || !fid_ref)
  {
    fprintf(stderr, "Error opening output files.\n");
    exit(EXIT_FAILURE);
  }

  // Write initial time tag
  uint64_t time_tag = 0x86ULL;
  fwrite(&time_tag, sizeof(uint64_t), 1, fid_out);

  long total_bytes = get_file_size(crystalID_file);
  if (total_bytes <= 0 || total_bytes % sizeof(uint32_t) != 0)
  {
    fprintf(stderr, "Invalid crystalID file size.\n");
    exit(EXIT_FAILURE);
  }
  long total_events = total_bytes / sizeof(uint32_t);
  long num_chunks = (total_events + CHUNK_EVENTS - 1) / CHUNK_EVENTS;
  long event_counter = 0;

  uint64_t last_cycle_time = 0;

  // Allocate buffers
  Event *events = malloc(CHUNK_EVENTS * sizeof(Event));
  if (!events)
  {
    fprintf(stderr, "Memory allocation failed.\n");
    exit(EXIT_FAILURE);
  }

  for (long ch = 0; ch < num_chunks; ch++)
  {
    long ne = (event_counter + CHUNK_EVENTS <= total_events) ? CHUNK_EVENTS : (total_events - event_counter);

    // Read data
    for (long i = 0; i < ne; i++)
    {
      fread(&events[i].crID, sizeof(uint32_t), 1, fid_cr);
      fread(&events[i].en, sizeof(float), 1, fid_en);
      fread(&events[i].evID, sizeof(uint32_t), 1, fid_ev);
      fread(&events[i].soID, sizeof(uint32_t), 1, fid_so);
      fread(&events[i].cpID, sizeof(uint32_t), 1, fid_cp);

      double tm_val;
      fread(&tm_val, sizeof(double), 1, fid_tm);
      events[i].tmBits = (uint64_t)llround(tm_val / CFD_RESOLUTION);
    }

    // Sort by tmBits
    qsort(events, ne, sizeof(Event), compare_events);

    // Prepare output buffers
    int maxw = ne + ne; // worst case: every event triggers a cycle change
    uint64_t *out_buf = malloc(maxw * sizeof(uint64_t));
    if (!out_buf)
    {
      fprintf(stderr, "Out buffer allocation failed.\n");
      exit(EXIT_FAILURE);
    }

    uint32_t *ref_buf = malloc(3 * ne * sizeof(uint32_t));
    if (!ref_buf)
    {
      fprintf(stderr, "Ref buffer allocation failed.\n");
      exit(EXIT_FAILURE);
    }

    int idx_out = 0;
    for (long i = 0; i < ne; i++)
    {
      uint64_t dtBits = events[i].tmBits - last_cycle_time;
      if ((int64_t)dtBits < 0 || dtBits >= TIME_TAG_INTERVAL_CFD_RESOLUTION)
      {
        last_cycle_time = events[i].tmBits;
        out_buf[idx_out++] = 0x86ULL;
        dtBits = 0;
      }

      double crystal_d = (double)events[i].crID;
      uint32_t trID = (uint32_t)fmod(crystal_d, NB_CRYSTALS_PER_RING_L);
      uint32_t axID = (uint32_t)floor(crystal_d / NB_CRYSTALS_PER_RING_L);

      int64_t eBit_raw = llround((events[i].en - ENERGY_OFFSET) / ENERGY_RESOLUTION);
      uint32_t eBits = (eBit_raw < 0) ? 0 : (eBit_raw >= MAX_EN ? MAX_EN - 1 : (uint32_t)eBit_raw);

      // 根据新的结构体布局重新编码 64-bit word
      // 0-7: TAG_SINGLES (0x82ULL)
      // 8-18: trID (11位)
      // 19-28: axID (10位)
      // 29-37: eBits (9位)
      // 38-63: dtBits (26位)
      uint64_t w = 0x82ULL |
                   ((uint64_t)trID << 8) |
                   ((uint64_t)axID << 19) |
                   ((uint64_t)eBits << 29) |
                   (dtBits << 38);

      out_buf[idx_out++] = w;

      // Fill reference buffer (3 uint32 per event)
      ref_buf[3 * i] = events[i].evID;
      ref_buf[3 * i + 1] = events[i].soID;
      ref_buf[3 * i + 2] = events[i].cpID;
    }

    // Write to files
    fwrite(out_buf, sizeof(uint64_t), idx_out, fid_out);
    fwrite(ref_buf, sizeof(uint32_t), 3 * ne, fid_ref);

    free(out_buf);
    free(ref_buf);

    event_counter += ne;
    printf("Chunk %ld/%ld encoded, events %ld-%ld.\n",
           ch + 1, num_chunks, event_counter - ne + 1, event_counter);
  }

  // Final time tag
  fwrite(&time_tag, sizeof(uint64_t), 1, fid_out);

  // Close files
  fclose(fid_cr);
  fclose(fid_tm);
  fclose(fid_en);
  fclose(fid_ev);
  fclose(fid_so);
  fclose(fid_cp);
  fclose(fid_out);
  fclose(fid_ref);
  free(events);

  printf("Encoding complete: %s and %s\n", output_dat, refere_dat);
}

int main(int argc, char *argv[])
{
  if (argc != 4)
  {
    fprintf(stderr, "Usage: %s <input_path> <output_path> <base_filename>\n", argv[0]);
    fprintf(stderr, "  Example: %s /data/input/ /data/output/ mydata\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  const char *input_path = argv[1];
  const char *output_path = argv[2];
  const char *base_name = argv[3];

  char crystalID_file[1024], time_file[1024], energy_file[1024];
  char eventID_file[1024], sourceID_file[1024], comptonP_file[1024];
  char output_dat[1024], refere_dat[1024];

  const char *crystalID_suff = ".crystalID";
  const char *time_suff = ".time";
  const char *energy_suff = ".energy";
  const char *eventID_suff = ".eventID";
  const char *sourceID_suff = ".sourceID";
  const char *comptonP_suff = ".comptonPhantom";
  const char *output_suff = ".dat";
  const char *refere_suff = ".reference";

  snprintf(crystalID_file, sizeof(crystalID_file), "%s%s%s", input_path, base_name, crystalID_suff);
  snprintf(time_file, sizeof(time_file), "%s%s%s", input_path, base_name, time_suff);
  snprintf(energy_file, sizeof(energy_file), "%s%s%s", input_path, base_name, energy_suff);
  snprintf(eventID_file, sizeof(eventID_file), "%s%s%s", input_path, base_name, eventID_suff);
  snprintf(sourceID_file, sizeof(sourceID_file), "%s%s%s", input_path, base_name, sourceID_suff);
  snprintf(comptonP_file, sizeof(comptonP_file), "%s%s%s", input_path, base_name, comptonP_suff);
  snprintf(output_dat, sizeof(output_dat), "%s%s%s", output_path, base_name, output_suff);
  snprintf(refere_dat, sizeof(refere_dat), "%s%s%s", output_path, base_name, refere_suff);

  encode_singles(crystalID_file, time_file, energy_file,
                 eventID_file, sourceID_file, comptonP_file,
                 output_dat, refere_dat);

  return 0;
}
