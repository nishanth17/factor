#ifndef FACTOR_ECM_NEON_SCREEN_H
#define FACTOR_ECM_NEON_SCREEN_H

typedef struct b4n_context b4n_context;
typedef struct { const unsigned char *data; unsigned size; } b4n_value;
typedef struct {
    double wall, cpu;
    unsigned long long checksum, pairs;
    int qos_set, qos_get, relative_priority;
    unsigned requested_qos, observed_qos;
    unsigned long long thread_id;
    int thread_id_status;
} b4n_observation;
typedef struct { unsigned char bytes[2][96]; unsigned size[2]; } b4n_output;

unsigned b4n_abi(void);
void b4n_qos(b4n_observation *out);
const char *b4n_gmp_version(void);
b4n_context *b4n_create(unsigned width, b4n_value modulus);
void b4n_destroy(b4n_context *ctx);
int b4n_set(b4n_context *ctx, unsigned index, const b4n_value values[4],
             unsigned active);
/* arm:0=mpn,1=scalar26,2=NEON26,3=production-mpz,
 * 4=hybrid-scalar26,5=hybrid-NEON26,6=hybrid-layout-only.
 * operation:0=multiply,1=square. No executor integration.
 */
int b4n_once(b4n_context *ctx, unsigned arm, unsigned operation, unsigned index);
int b4n_export(b4n_context *ctx, b4n_output *output);
/* Timing is available only to the separately admitted Stage A child. */
int b4n_batch(b4n_context *ctx, unsigned arm, unsigned operation,
               unsigned repeats, unsigned offset, b4n_observation *out);
int b4n_validate(b4n_context *ctx, const b4n_observation *out);
#endif
