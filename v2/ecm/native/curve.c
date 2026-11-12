/* Included by packet.c to share the independently checked point primitives.
 * Canonical job state and speculative arithmetic are deliberately distinct.
 * Checkpoints export only canonical state; interrupted compound actions restart
 * from their saved input, retaining their complete already-incurred charges.
 */
#define CURVE_TABLE_LIMIT 65536u
#define CURVE_BLOCK_LIMIT 4096u
enum { CS_SETUP, CS_ONE, CS_REPLAY, CS_TWO_SETUP, CS_BABY, CS_GIANT,
       CS_TWO, CS_TERMS };
enum { CO_IDLE, CO_FAST, CO_FAST_CHECK, CO_RECOVER, CO_STRICT,
       CO_STRICT_CHECK, CO_REPLAY, CO_REPLAY_CHECK, CO_GIANT_A,
       CO_GIANT_B, CO_SMALL, CO_GROUP };

struct b4_prime_block {
    uint64_t lo, hi, bound;
    unsigned chain, refs;
    size_t count;
    b4_curve_unit *units;
    uint64_t b1, b2, *centers;
    unsigned *offsets, prepared;
};

struct b4_curve_state {
    mpz_t value[2], replay[2], giant[2], previous[2], product, factor;
    mpz_t temporary[2], recover_value[2], unit_saved[2], strict_saved[2];
    mpz_t temp_product, divisor, terms[256];
    mpz_t *bx, *bz;
    size_t table_capacity;
    b4_curve_unit powers[256];
    b4_record records[256];
    b4_prime_block *block, *blocks[CURVE_BLOCK_LIMIT];
    size_t block_count;
    b4_curve_info info;
    uint64_t sigma, remaining, meter;
    size_t recover_index, strict_opcode, fast_completed, group_end, pending;
    unsigned active, chunk, batch, segment, operation, strict_phase;
    unsigned strict_result, unit_mode, atomic_actions, phase_clocks;
    long strict_bit;
    const b4_record *strict_record;
    uint64_t strict_scalar;
    atomic_uint cancelled;
    unsigned preparation_failure, cancel_after_actions;
};

unsigned b4_whole_ecm(void) { return 1; }

unsigned b4_coalesced_unit_yields(void) {
#ifdef B4_COALESCE_UNIT_YIELDS
    return 1;
#else
    return 0;
#endif
}

static unsigned bits64(uint64_t value) {
    unsigned bits = 0;
    do { ++bits; value >>= 1; } while (value);
    return bits;
}

static uint64_t sqrt64(uint64_t value) {
    uint64_t lo = 0, hi = (uint64_t)1 << 32;
    while (lo + 1 < hi) {
        uint64_t mid = lo + (hi - lo) / 2;
        if (mid <= value / mid) lo = mid; else hi = mid;
    }
    return lo;
}

static int record_valid(const b4_record *r, uint64_t scalar) {
    if (!r->code || !r->masks || r->code_size % 4 || r->code_size > 2048 ||
        r->masks_size != r->code_size / 4 + 1 || !r->slots ||
        r->slots > POINTS || r->output >= r->slots || r->masks[0] > 3 ||
        !scalar || scalar > UINT32_MAX) return 0;
    uint64_t values[POINTS] = {1};
    for (size_t j = 0; j < r->code_size; j += 4) {
        const unsigned char *op = r->code + j;
        if (op[0] >= r->slots || op[1] >= r->slots || op[2] >= r->slots ||
            !values[op[1]] || !values[op[2]] || r->masks[j / 4 + 1] > 3)
            return 0;
        uint64_t a = values[op[1]], b = values[op[2]], v;
        if (op[3] == 255) {
            if (op[1] != op[2]) return 0;
            v = 2 * a;
        } else {
            if (op[3] >= r->slots || !values[op[3]]) return 0;
            uint64_t difference = a > b ? a - b : b - a;
            if (values[op[3]] == difference) v = a + b;
            else if (values[op[3]] == a + b) v = difference;
            else return 0;
        }
        if (!v || v > 2 * scalar + 2) return 0;
        values[op[0]] = v;
    }
    return values[r->output] == scalar;
}

static int record_copy(b4_record *dest, const b4_record *src) {
    unsigned char *code = malloc(src->code_size ? src->code_size : 1);
    unsigned char *masks = malloc(src->masks_size);
    if (!code || !masks) { free(code); free(masks); return -1; }
    memcpy(code, src->code, src->code_size);
    memcpy(masks, src->masks, src->masks_size);
    *dest = *src; dest->code = code; dest->masks = masks;
    return 0;
}

/* This is verify_frontier's forward-ancestor certificate, independently
 * checked before native-only catalogue use. At most512 instructions introduce
 * 1026 coordinate bits; fixed scratch is bounded independently of the input. */
static int frontier_valid(const b4_record *r, uint64_t scalar) {
    if (!record_valid(r, scalar)) return 0;
    uint64_t states[POINTS][2][17] = {{{0}}}, covered[17] = {0};
    states[0][0][0] = 1; states[0][1][0] = 2;
    covered[0] = r->masks[0];
    size_t steps = r->code_size / 4;
    for (size_t index = 1; index <= steps; ++index) {
        const unsigned char *op = r->code + 4 * (index - 1);
        uint64_t x[17] = {0}, z[17] = {0};
        x[(2 * index) / 64] = (uint64_t)1 << ((2 * index) % 64);
        z[(2 * index + 1) / 64] = (uint64_t)1 << ((2 * index + 1) % 64);
        for (unsigned word = 0; word < 17; ++word) {
            if (op[3] == 255)
                z[word] |= states[op[1]][0][word] | states[op[1]][1][word];
            else {
                x[word] |= states[op[3]][1][word];
                z[word] |= states[op[3]][0][word];
            }
            if (r->masks[index] & 1) covered[word] |= x[word];
            if (r->masks[index] & 2) covered[word] |= z[word];
        }
        memcpy(states[op[0]][0], x, sizeof(x));
        memcpy(states[op[0]][1], z, sizeof(z));
    }
    unsigned bits = (unsigned)(2 * (steps + 1));
    for (unsigned word = 0; word < 17; ++word) {
        uint64_t actual = covered[word] | states[r->output][0][word] |
                          states[r->output][1][word];
        unsigned left = bits > 64 * word ? bits - 64 * word : 0;
        uint64_t expected = left >= 64 ? UINT64_MAX :
                            left ? ((uint64_t)1 << left) - 1 : 0;
        if (actual != expected) return 0;
    }
    return 1;
}

struct b4_catalogue {
    size_t count;
    uint64_t *scalars;
    b4_record *records;
};

static uint64_t little_get(const unsigned char *p, unsigned width) {
    uint64_t value = 0;
    for (unsigned i = 0; i < width; ++i) value |= (uint64_t)p[i] << (8 * i);
    return value;
}

void b4_catalogue_destroy(b4_catalogue *catalogue) {
    if (!catalogue) return;
    for (size_t i = 0; i < catalogue->count; ++i) {
        free((void *)catalogue->records[i].code);
        free((void *)catalogue->records[i].masks);
    }
    free(catalogue->records); free(catalogue->scalars); free(catalogue);
}

b4_catalogue *b4_catalogue_create(const unsigned char *data, size_t size) {
    if (!data || size < 4 || size > 2 * 1024 * 1024) return NULL;
    size_t count = (size_t)little_get(data, 4), offset = 4;
    if (count > 1024) return NULL;
    b4_catalogue *catalogue = calloc(1, sizeof(*catalogue));
    if (!catalogue) return NULL;
    catalogue->scalars = calloc(count ? count : 1, sizeof(uint64_t));
    catalogue->records = calloc(count ? count : 1, sizeof(b4_record));
    if (!catalogue->scalars || !catalogue->records) goto invalid_catalogue;
    uint64_t previous = 0;
    for (size_t i = 0; i < count; ++i) {
        if (size - offset < 14) goto invalid_catalogue;
        const unsigned char *p = data + offset;
        uint64_t scalar = little_get(p, 8);
        size_t codes = (size_t)little_get(p + 8, 2);
        size_t masks = (size_t)little_get(p + 10, 2);
        offset += 14;
        if (scalar <= previous || codes > size - offset ||
            masks > size - offset - codes) goto invalid_catalogue;
        b4_record record = {data + offset, data + offset + codes,
                            codes, masks, p[12], p[13]};
        if (!frontier_valid(&record, scalar) ||
            record_copy(&catalogue->records[i], &record)) goto invalid_catalogue;
        catalogue->scalars[i] = scalar;
        ++catalogue->count; previous = scalar; offset += codes + masks;
    }
    if (offset != size) goto invalid_catalogue;
    return catalogue;
invalid_catalogue:
    b4_catalogue_destroy(catalogue);
    return NULL;
}

static const b4_record *catalogue_find(const b4_catalogue *catalogue,
                                       uint64_t scalar) {
    size_t lo = 0, hi = catalogue->count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (catalogue->scalars[mid] < scalar) lo = mid + 1; else hi = mid;
    }
    return lo < catalogue->count && catalogue->scalars[lo] == scalar ?
           &catalogue->records[lo] : NULL;
}

/* The same finite segmented Eratosthenes contract as SieveContext. This is a
 * bulk implementation of existing work, not a different prime algorithm.
 * Wider endpoints retain Python's exact implementation. The root is <=65535,
 * block width <=65536, and every retained descriptor is checked before use. */
b4_prime_block *b4_block_generate(b4_catalogue *catalogue,
                                  unsigned long long lo, unsigned long long hi,
                                  unsigned long long bound, unsigned chain) {
    if (!catalogue || lo < 2 || hi <= lo || hi > ((uint64_t)1 << 32) ||
        hi - lo > 65536 || chain > 1 || (chain && (!bound || bound > 2000)) ||
        (bound && hi - 1 > bound)) return NULL;
    size_t width = (size_t)(hi - lo), root = (size_t)sqrt64(hi - 1);
    unsigned char *base = calloc(root + 1, 1), *composite = calloc(width, 1);
    b4_curve_unit *units = calloc(width / 2 + 1, sizeof(*units));
    if (!base || !composite || !units) {
        free(base); free(composite); free(units); return NULL;
    }
    for (uint64_t prime = 2; prime <= root; ++prime) {
        if (base[prime]) continue;
        if (prime <= root / prime)
            for (uint64_t multiple = prime * prime; multiple <= root;
                 multiple += prime) base[multiple] = 1;
        uint64_t first = ((lo + prime - 1) / prime) * prime;
        if (first < prime * prime) first = prime * prime;
        for (uint64_t multiple = first; multiple < hi; multiple += prime)
            composite[multiple - lo] = 1;
    }
    size_t count = 0;
    unsigned valid = 1;
    for (uint64_t prime = lo; prime < hi; ++prime) {
        if (composite[prime - lo]) continue;
        b4_curve_unit *u = &units[count++];
        u->prime = u->power = prime;
        if (bound) while (u->power <= bound / prime) u->power *= prime;
        if (chain) {
            const b4_record *fast = catalogue_find(catalogue, u->power);
            const b4_record *unit = catalogue_find(catalogue, prime);
            if (!fast || !unit) { valid = 0; break; }
            u->fast = *fast; u->unit = *unit;
            unsigned exponent = 0;
            for (uint64_t power = u->power; power > 1; power /= prime) ++exponent;
            /* Match ChainPlan's funded fast/strict/repeated-unit work. */
            u->reserve = 3 * (fast->code_size / 4) + 3 +
                5 * (fast->code_size / 4) + 10 * bits64(u->power) + 8 +
                exponent * (5 * (unit->code_size / 4) + 10 * bits64(prime) + 8);
        }
    }
    b4_prime_block *block = valid ?
        b4_block_create(lo, hi, bound, chain, units, count) : NULL;
    free(base); free(composite); free(units);
    return block;
}

size_t b4_block_count(const b4_prime_block *block) {
    return block ? block->count : 0;
}

size_t b4_block_bytes(const b4_prime_block *block) {
    if (!block) return 0;
    size_t bytes = sizeof(*block) + (block->count ? block->count : 1) * sizeof(b4_curve_unit);
    if (block->centers) bytes += block->count * sizeof(uint64_t);
    if (block->offsets) bytes += block->count * sizeof(unsigned);
    if (block->chain)
        for (size_t i = 0; i < block->count; ++i)
            bytes += block->units[i].fast.code_size + block->units[i].fast.masks_size +
                     block->units[i].unit.code_size + block->units[i].unit.masks_size;
    return bytes;
}

int b4_block_encode(const b4_prime_block *block, unsigned char *data, size_t size) {
    if (!block || size != 16 * block->count || (size && !data)) return -1;
    for (size_t i = 0; i < block->count; ++i) {
        uint64_t values[2] = {block->units[i].prime, block->units[i].power};
        for (unsigned field = 0; field < 2; ++field)
            for (unsigned byte = 0; byte < 8; ++byte)
                data[16 * i + 8 * field + byte] =
                    (unsigned char)(values[field] >> (8 * byte));
    }
    return 0;
}

b4_prime_block *b4_block_decode(b4_catalogue *catalogue,
                                unsigned long long lo, unsigned long long hi,
                                unsigned long long bound, unsigned chain,
                                const unsigned char *data, size_t size) {
    if ((size && !data) || size % 16 || size > 65536 * 16) return NULL;
    b4_prime_block *block = b4_block_generate(catalogue, lo, hi, bound, chain);
    if (!block) return NULL;
    /* Hashes identify artifacts; this independently proves exact coverage and
     * powers, including no missing/extra prime. Programs come only from the
     * bound, certified immutable catalogue, never unchecked cached opcodes. */
    unsigned valid = size == 16 * block->count;
    for (size_t i = 0; valid && i < block->count; ++i)
        valid = little_get(data + 16 * i, 8) == block->units[i].prime &&
                little_get(data + 16 * i + 8, 8) == block->units[i].power;
    if (!valid) { b4_block_destroy(block); return NULL; }
    return block;
}

int b4_block_destroy(b4_prime_block *b) {
    if (!b) return 0;
    if (b->refs) return -1;
    for (size_t i = 0; i < b->count; ++i) {
        free((void *)b->units[i].fast.code); free((void *)b->units[i].fast.masks);
        free((void *)b->units[i].unit.code); free((void *)b->units[i].unit.masks);
    }
    free(b->units); free(b->centers); free(b->offsets); free(b);
    return 0;
}

b4_prime_block *b4_block_create(unsigned long long lo, unsigned long long hi,
                               unsigned long long bound, unsigned chain,
                               const b4_curve_unit *units, size_t count) {
    if (lo < 2 || hi <= lo || count > 65536 || (count && !units) || chain > 1 ||
        (bound && (bound < 2 || hi - 1 > bound)) || (!bound && chain)) return NULL;
    uint64_t previous = lo - 1;
    for (size_t i = 0; i < count; ++i) {
        const b4_curve_unit *u = &units[i];
        if (u->prime <= previous || u->prime < lo || u->prime >= hi ||
            (u->prime != 2 && !(u->prime & 1))) return NULL;
        uint64_t power = u->prime;
        if (bound) while (power <= bound / u->prime) power *= u->prime;
        if (u->power != power || (chain &&
            (!record_valid(&u->fast, power) || !record_valid(&u->unit, u->prime))))
            return NULL;
        previous = u->prime;
    }
    b4_prime_block *b = calloc(1, sizeof(*b));
    if (!b) return NULL;
    b->units = calloc(count ? count : 1, sizeof(*b->units));
    if (!b->units) { free(b); return NULL; }
    b->lo = lo; b->hi = hi; b->bound = bound; b->chain = chain; b->count = count;
    for (size_t i = 0; i < count; ++i) {
        b->units[i].prime = units[i].prime; b->units[i].power = units[i].power;
        b->units[i].reserve = units[i].reserve;
        if (chain && (record_copy(&b->units[i].fast, &units[i].fast) ||
                      record_copy(&b->units[i].unit, &units[i].unit))) {
            b4_block_destroy(b); return NULL;
        }
    }
    return b;
}

int b4_block_prepare(b4_prime_block *b, unsigned long long b1,
                      unsigned long long b2) {
    if (!b || b->bound || b1 < 2 || b2 < b1 || b2 == UINT64_MAX ||
        b->lo < b1 + 1 || b->hi > b2 + 1) return -1;
    if (b->prepared) return b->b1 == b1 && b->b2 == b2 ? 0 : -1;
    uint64_t origin = b1 & 1 ? b1 : b1 - 1, distance = sqrt64(b2);
    if (distance > (origin - 1) / 2) distance = (origin - 1) / 2;
    if (distance > CURVE_TABLE_LIMIT) return -2;
    if (distance >= 2 && b->count) {
        uint64_t *centers = calloc(b->count, sizeof(*centers));
        unsigned *offsets = calloc(b->count, sizeof(*offsets));
        if (!centers || !offsets) { free(centers); free(offsets); return -5; }
        for (size_t i = 0; i < b->count; ++i) {
            uint64_t p = b->units[i].prime;
            uint64_t group = (p - origin - 1) / (2 * distance);
            centers[i] = origin + group * 2 * distance;
            offsets[i] = (unsigned)((p - centers[i]) / 2);
            if (!offsets[i] || offsets[i] > distance) { free(centers); free(offsets); return -1; }
        }
        b->centers = centers; b->offsets = offsets;
    }
    b->b1 = b1; b->b2 = b2; b->prepared = 1;
    return 0;
}

static void curve_abort(b4_context *c) {
    if (!c || !c->curve) return;
    b4_curve_state *w = c->curve;
    for (size_t i = 0; i < w->block_count; ++i) --w->blocks[i]->refs;
    w->block_count = 0; w->block = NULL; w->active = 0;
}

static void curve_free(b4_context *c) {
    if (!c || !c->curve) return;
    curve_abort(c);
    b4_curve_state *w = c->curve;
    for (unsigned i = 0; i < 2; ++i)
        mpz_clears(w->value[i], w->replay[i], w->giant[i], w->previous[i],
                   w->temporary[i], w->recover_value[i], w->unit_saved[i],
                   w->strict_saved[i], NULL);
    for (unsigned i = 0; i < 256; ++i) mpz_clear(w->terms[i]);
    for (size_t i = 0; i < w->table_capacity; ++i) {
        mpz_clear(w->bx[i]); mpz_clear(w->bz[i]);
    }
    free(w->bx); free(w->bz);
    mpz_clears(w->product, w->factor, w->temp_product, w->divisor, NULL);
    free(w); c->curve = NULL;
}

static int curve_active(const b4_context *c) {
    return c && c->curve && c->curve->active;
}

static int curve_allocate(b4_context *c) {
    if (c->curve) return 0;
    b4_curve_state *w = calloc(1, sizeof(*w));
    if (!w) return -5;
    unsigned bits = (unsigned)mpz_sizeinbase(c->n, 2);
    for (unsigned i = 0; i < 2; ++i) {
        mpz_init2(w->value[i], bits); mpz_init2(w->replay[i], bits);
        mpz_init2(w->giant[i], bits); mpz_init2(w->previous[i], bits);
        mpz_init2(w->temporary[i], bits); mpz_init2(w->recover_value[i], bits);
        mpz_init2(w->unit_saved[i], bits); mpz_init2(w->strict_saved[i], bits);
    }
    for (unsigned i = 0; i < 256; ++i) mpz_init2(w->terms[i], bits);
    mpz_init2(w->product, bits); mpz_init2(w->factor, bits);
    mpz_init2(w->temp_product, 3 * bits + 64); mpz_init2(w->divisor, bits);
    c->curve = w;
    atomic_init(&w->cancelled, 0);
    return 0;
}

static int curve_table(b4_context *c, size_t count) {
    b4_curve_state *w = c->curve;
    if (w->preparation_failure) return -5;
    if (count <= w->table_capacity) return 0;
    if (count > CURVE_TABLE_LIMIT + 1) return -2;
    mpz_t *x = calloc(count, sizeof(mpz_t)), *z = calloc(count, sizeof(mpz_t));
    if (!x || !z) { free(x); free(z); return -5; }
    unsigned bits = (unsigned)mpz_sizeinbase(c->n, 2);
    for (size_t i = 0; i < count; ++i) { mpz_init2(x[i], bits); mpz_init2(z[i], bits); }
    for (size_t i = 0; i < w->table_capacity; ++i) {
        mpz_set(x[i], w->bx[i]); mpz_set(z[i], w->bz[i]);
        mpz_clear(w->bx[i]); mpz_clear(w->bz[i]);
    }
    free(w->bx); free(w->bz); w->bx = x; w->bz = z; w->table_capacity = count;
    return 0;
}

int b4_curve_begin(b4_context *c, b4_value n, b4_value sigma,
                    unsigned long long b1, unsigned long long b2,
                    unsigned chunk, unsigned batch, unsigned segment) {
    if (c) b4_abort(c);
    if (!c || !canonical(n, 128) || !canonical(sigma, 8) || b1 < 2 || b2 < b1 ||
        b2 == UINT64_MAX || !chunk || chunk > 256 || !batch || batch > 256 ||
        !segment || segment > 32768) return -1;
    uint64_t distance = sqrt64(b2), center = b1 & 1 ? b1 : b1 - 1;
    if (distance > (center - 1) / 2) distance = (center - 1) / 2;
    if (distance > CURVE_TABLE_LIMIT) return -2;
    curve_abort(c);
    import_value(c->n, n); import_value(c->scalar, sigma);
    if (mpz_cmp_ui(c->n, 3) < 0 || mpz_even_p(c->n) || mpz_cmp_ui(c->scalar, 6) < 0)
        return -1;
    unsigned bits = (unsigned)mpz_sizeinbase(c->n, 2);
    if (bits > c->capacity_bits) {
        for (unsigned i = 0; i < POINTS; ++i) {
            mpz_realloc2(c->x[i], bits); mpz_realloc2(c->z[i], bits);
        }
        for (unsigned i = 0; i < 7; ++i) mpz_realloc2(c->scratch[i], 4 * bits + 64);
        mpz_realloc2(c->a24, bits); mpz_realloc2(c->guard, 3 * bits + 64);
        mpz_realloc2(c->wide_product, 4 * bits + 64);
        mpz_realloc2(c->scalar, 16384); c->capacity_bits = bits;
    }
    int status = curve_allocate(c);
    if (status) return status;
    b4_curve_state *w = c->curve;
    memset(&w->info, 0, sizeof(w->info));
    w->info.b1 = b1; w->info.b2 = b2; w->info.phase = CS_SETUP;
    w->info.next = w->info.left = 2; w->info.hi = b1 + 1;
    w->sigma = mpz_get_ui(c->scalar); w->active = 1;
    w->chunk = chunk; w->batch = batch; w->segment = segment;
    w->operation = CO_IDLE; w->pending = 0; w->strict_phase = 0;
    atomic_store(&w->cancelled, 0);
    w->preparation_failure = w->cancel_after_actions = 0;
    mpz_set_ui(w->factor, 0); mpz_set_ui(w->product, 1);
    c->mode = c->done = 0;
    /* Legacy packet rebinding must not trust stale cached modulus/a24 bytes. */
    c->modulus_size = c->a24_size = 0;
    return 0;
}

int b4_curve_options(b4_context *c, unsigned atomic_actions, unsigned clocks) {
    if (!curve_active(c) || atomic_actions > 1 || clocks > 1) return -1;
    c->curve->atomic_actions = atomic_actions; c->curve->phase_clocks = clocks;
    return 0;
}

int b4_curve_faults(b4_context *c, unsigned preparation_failure,
                     unsigned cancel_after_actions) {
    if (!curve_active(c) || preparation_failure > 1) return -1;
    c->curve->preparation_failure = preparation_failure;
    c->curve->cancel_after_actions = cancel_after_actions;
    return 0;
}

void b4_curve_cancel(b4_context *c) {
    if (c && c->curve) atomic_store(&c->curve->cancelled, 1);
}

int b4_curve_feed(b4_context *c, b4_prime_block *b) {
    if (!curve_active(c) || !b) return -1;
    b4_curve_state *w = c->curve;
    uint64_t end = w->info.next + 2 * (uint64_t)w->segment;
    if (end < w->info.next || end > w->info.hi) end = w->info.hi;
    uint64_t bound = w->info.phase == CS_ONE ? w->info.b1 : 0;
    if (!w->info.need || b->lo != w->info.next || b->hi != end || b->bound != bound)
        return -1;
    if (!bound && (!b->prepared || b->b1 != w->info.b1 || b->b2 != w->info.b2)) return -1;
    size_t i;
    for (i = 0; i < w->block_count; ++i) if (w->blocks[i] == b) break;
    if (i == w->block_count) {
        if (i >= CURVE_BLOCK_LIMIT) return -2;
        w->blocks[w->block_count++] = b; ++b->refs;
    }
    w->block = b; w->info.left = b->lo; w->info.next = b->hi;
    w->info.index = 0; w->info.need = 0;
    return 0;
}

static void pair_set(mpz_t *dest, mpz_t *src) {
    mpz_set(dest[0], src[0]); mpz_set(dest[1], src[1]);
}
static void load_pair(b4_context *c, unsigned slot, mpz_t *point) {
    mpz_set(c->x[slot], point[0]); mpz_set(c->z[slot], point[1]);
}
static void save_pair(b4_context *c, mpz_t *point, unsigned slot) {
    mpz_set(point[0], c->x[slot]); mpz_set(point[1], c->z[slot]);
}
static void math_begin(b4_context *c, mpz_t *point) {
    load_pair(c, 0, point); load_pair(c, 1, point);
    c->mode = 1; c->done = 0;
    if (!mpz_sgn(c->scalar) || !mpz_sgn(c->z[0])) {
        mpz_set_ui(c->x[1], 1); mpz_set_ui(c->z[1], 0); c->done = 1;
    } else if (mpz_cmp_ui(c->scalar, 1) == 0) c->done = 1;
    else { doubling(c, 2, 0); c->bit = (long)mpz_sizeinbase(c->scalar, 2) - 2; }
}
static void math_ui(b4_context *c, mpz_t *point, uint64_t scalar) {
    mpz_set_ui(c->scalar, scalar); math_begin(c, point);
}
static void curve_finish(b4_context *c, const mpz_t factor) {
    b4_curve_state *w = c->curve;
    w->info.done = 1;
    if (factor) mpz_set(w->factor, factor); else mpz_set_ui(w->factor, 0);
    w->info.boundary = 1;
}
static int proper(b4_context *c, const mpz_t factor) {
    return mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, c->n) < 0;
}

/* Z first, then X, matching _check_point and Executor's explicit X guard.
 * 0=certified, 1=proper factor, 2=exceptional (checked-ladder retry required). */
static int strict_check(b4_context *c, unsigned slot) {
    b4_curve_state *w = c->curve;
    mpz_gcd(w->divisor, c->z[slot], c->n);
    if (proper(c, w->divisor)) { curve_finish(c, w->divisor); return 1; }
    if (mpz_cmp(w->divisor, c->n) == 0) {
        mpz_gcd(w->divisor, c->x[slot], c->n);
        if (proper(c, w->divisor)) { curve_finish(c, w->divisor); return 1; }
        if (mpz_cmp(w->divisor, c->n) == 0) return 2;
    }
    mpz_gcd(w->divisor, c->x[slot], c->n);
    if (proper(c, w->divisor)) { curve_finish(c, w->divisor); return 1; }
    return 0;
}

static void strict_begin(b4_context *c, const b4_record *r, uint64_t scalar) {
    b4_curve_state *w = c->curve;
    w->strict_record = r; w->strict_scalar = scalar;
    pair_set(w->strict_saved, w->recover_value);
    load_pair(c, 0, w->recover_value);
    w->strict_opcode = 0; w->strict_phase = 1; w->strict_result = 0;
    w->operation = CO_STRICT;
}

/* One strict point operation per coroutine action, retaining register state.
 * Its canonical chunk input remains unchanged until the whole funded compound
 * action completes, including the exceptional-difference ladder and unit replay. */
static void strict_action(b4_context *c) {
    b4_curve_state *w = c->curve;
    const b4_record *r = w->strict_record;
    if (w->strict_phase == 1) {
        int checked = strict_check(c, 0);
        if (checked == 1) return;
        if (!checked && mpz_sgn(c->x[0]) && mpz_sgn(c->z[0])) w->strict_phase = 2;
        else w->strict_phase = 3;
    } else if (w->strict_phase == 2) {
        if (w->strict_opcode == r->code_size) {
            save_pair(c, w->recover_value, r->output); w->strict_result = 1;
        } else {
            const unsigned char *op = r->code + w->strict_opcode;
            if (op[3] != 255 && (!mpz_sgn(c->x[op[3]]) || !mpz_sgn(c->z[op[3]]))) {
                w->strict_phase = 3; return;
            }
            if (op[3] == 255) doubling(c, op[0], op[1]);
            else addition(c, op[0], op[1], op[2], op[3]);
            int checked = strict_check(c, op[0]);
            if (checked == 1) return;
            if (checked == 2) w->strict_phase = 3;
            else w->strict_opcode += 4;
        }
    } else if (w->strict_phase == 3) {
        load_pair(c, 0, w->strict_saved); load_pair(c, 1, w->strict_saved);
        if (!mpz_sgn(c->z[0]) || !w->strict_scalar ||
            (!mpz_sgn(c->x[0]) && !(w->strict_scalar & 1))) {
            mpz_set_ui(w->recover_value[0], 1); mpz_set_ui(w->recover_value[1], 0);
            w->strict_result = 1;
        } else if (!mpz_sgn(c->x[0]) || w->strict_scalar == 1) {
            pair_set(w->recover_value, w->strict_saved); w->strict_result = 1;
        } else {
            doubling(c, 2, 0);
            int checked = strict_check(c, 2);
            if (checked == 1) return;
            if (checked == 2) w->strict_result = 2;
            else { w->strict_bit = (long)bits64(w->strict_scalar) - 2; w->strict_phase = 4; }
        }
    } else if (w->strict_phase == 4) {
        addition(c, 3, 1, 2, 0);
        int checked = strict_check(c, 3);
        if (checked == 1) return;
        if (checked == 2) w->strict_result = 2; else w->strict_phase = 5;
    } else {
        unsigned bit = (unsigned)((w->strict_scalar >> w->strict_bit) & 1);
        doubling(c, 4, bit ? 2 : 1);
        int checked = strict_check(c, 4);
        if (checked == 1) return;
        if (checked == 2) w->strict_result = 2;
        else {
            mpz_set(c->x[1], c->x[bit ? 3 : 4]); mpz_set(c->z[1], c->z[bit ? 3 : 4]);
            mpz_set(c->x[2], c->x[bit ? 4 : 3]); mpz_set(c->z[2], c->z[bit ? 4 : 3]);
            if (--w->strict_bit < 0) { save_pair(c, w->recover_value, 1); w->strict_result = 1; }
            else w->strict_phase = 4;
        }
    }
    if (w->strict_result) w->operation = CO_STRICT_CHECK;
}

static void chunk_replay(b4_context *c) {
    b4_curve_state *w = c->curve;
    if (w->info.chain) { ++w->info.chain_chunks; ++w->info.chain_replays; }
    w->info.phase = CS_REPLAY; w->info.replay_index = 0; w->info.replay_power = 1;
    pair_set(w->replay, w->value); w->operation = CO_IDLE; w->info.boundary = 1;
}

static void chunk_commit(b4_context *c, mpz_t *point, unsigned replayed) {
    b4_curve_state *w = c->curve;
    pair_set(w->value, point); w->info.powers_count = 0;
    if (w->info.chain) { ++w->info.chain_chunks; w->info.chain_replays += replayed; }
    w->operation = CO_IDLE; w->info.boundary = 1;
}

static int setup_action(b4_context *c) {
    b4_curve_state *w = c->curve;
    mpz_t *s = c->scratch;
    w->info.work += mpz_sizeinbase(c->n, 2);
    mpz_set_ui(s[0], w->sigma); mpz_mul(s[0], s[0], s[0]);
    mpz_sub_ui(s[0], s[0], 5); mpz_mod(s[0], s[0], c->n);
    mpz_set_ui(s[1], w->sigma); mpz_mul_ui(s[1], s[1], 4); mpz_mod(s[1], s[1], c->n);
    mpz_powm_ui(w->temporary[0], s[0], 3, c->n);
    mpz_powm_ui(w->temporary[1], s[1], 3, c->n);
    mpz_mul(s[2], w->temporary[0], s[1]); mpz_mul_ui(s[2], s[2], 16);
    mpz_mod(s[2], s[2], c->n); mpz_gcd(w->divisor, s[2], c->n);
    if (proper(c, w->divisor)) { curve_finish(c, w->divisor); return 0; }
    if (mpz_cmp(w->divisor, c->n) == 0) { curve_finish(c, NULL); return 0; }
    if (!mpz_invert(s[3], s[2], c->n)) return -1;
    mpz_sub(s[2], s[1], s[0]); mpz_powm_ui(s[2], s[2], 3, c->n);
    mpz_mul_ui(s[4], s[0], 3); mpz_add(s[4], s[4], s[1]);
    mpz_mul(s[2], s[2], s[4]); mpz_mul(s[2], s[2], s[3]);
    mpz_mod(c->a24, s[2], c->n);
    mpz_mul_ui(s[2], c->a24, 4); mpz_sub_ui(s[2], s[2], 2);
    mpz_mod(s[2], s[2], c->n); mpz_mul(s[2], s[2], s[2]); mpz_sub_ui(s[2], s[2], 4);
    mpz_gcd(w->divisor, s[2], c->n);
    if (proper(c, w->divisor)) { curve_finish(c, w->divisor); return 0; }
    if (mpz_cmp(w->divisor, c->n) == 0) { curve_finish(c, NULL); return 0; }
    pair_set(w->value, w->temporary); w->info.setup_done = 1;
    w->info.phase = CS_ONE; w->info.boundary = 1;
    return 0;
}

static int need_prime(b4_curve_state *w) {
    if (w->block && w->info.index < w->block->count) return 1;
    if (w->info.next < w->info.hi) { w->info.need = 1; return 2; }
    return 0;
}

static int curve_action(b4_context *c) {
    b4_curve_state *w = c->curve;
    b4_curve_info *i = &w->info;
    if (w->operation == CO_FAST || w->operation == CO_REPLAY ||
        w->operation == CO_GIANT_A || w->operation == CO_GIANT_B ||
        w->operation == CO_SMALL) {
        if (!c->done) arithmetic_action(c);
        if (!c->done) return 0;
        unsigned slot = c->mode == 1 ? 1 : 0;
        save_pair(c, w->temporary, slot);
        if (w->operation == CO_FAST) { w->fast_completed = c->completed; w->operation = CO_FAST_CHECK; }
        else if (w->operation == CO_REPLAY) w->operation = CO_REPLAY_CHECK;
        else if (w->operation == CO_GIANT_A) {
            pair_set(w->recover_value, w->temporary);
            math_ui(c, w->value, i->center - 2 * i->distance); w->operation = CO_GIANT_B;
        } else if (w->operation == CO_GIANT_B) {
            pair_set(w->giant, w->recover_value); pair_set(w->previous, w->temporary);
            i->phase = CS_TWO; w->operation = CO_IDLE; i->boundary = 1;
        } else {
            mpz_set(w->terms[i->terms_count], w->temporary[1]);
            modular_product(c, w->product, w->product, w->temporary[1]);
            ++i->terms_count; ++i->index; w->operation = CO_IDLE; i->boundary = 1;
        }
        return 0;
    }
    if (w->operation == CO_FAST_CHECK) {
        if (!i->chain) {
            mpz_gcd(w->divisor, w->temporary[1], c->n);
            if (proper(c, w->divisor)) curve_finish(c, w->divisor);
            else if (mpz_cmp(w->divisor, c->n) == 0) chunk_replay(c);
            else chunk_commit(c, w->temporary, 0);
        } else {
            mpz_mul(w->divisor, c->guard, w->temporary[0]);
            mpz_mul(w->divisor, w->divisor, w->temporary[1]); mpz_gcd(w->divisor, w->divisor, c->n);
            if (mpz_cmp_ui(w->divisor, 1) == 0) chunk_commit(c, w->temporary, 0);
            else { pair_set(w->recover_value, w->value); w->recover_index = 0; w->unit_mode = 0; w->operation = CO_RECOVER; }
        }
        return 0;
    }
    if (w->operation == CO_RECOVER) {
        if (w->recover_index >= w->fast_completed) {
            if (w->fast_completed == i->powers_count) chunk_commit(c, w->recover_value, 1);
            else chunk_replay(c);
        } else {
            const b4_curve_unit *u = &w->powers[w->recover_index];
            if (!w->unit_mode) pair_set(w->unit_saved, w->recover_value);
            strict_begin(c, w->unit_mode ? &u->unit : &u->fast, w->unit_mode ? u->prime : u->power);
        }
        return 0;
    }
    if (w->operation == CO_STRICT) { strict_action(c); return 0; }
    if (w->operation == CO_STRICT_CHECK) {
        const b4_curve_unit *u = &w->powers[w->recover_index];
        mpz_gcd(w->divisor, w->recover_value[1], c->n);
        if (!w->unit_mode) {
            if (w->strict_result == 1 && mpz_cmp(w->divisor, c->n)) ++w->recover_index;
            else { pair_set(w->recover_value, w->unit_saved); w->remaining = u->power; w->unit_mode = 1; }
        } else {
            if (w->strict_result != 1 || mpz_cmp(w->divisor, c->n) == 0) { chunk_replay(c); return 0; }
            w->remaining /= u->prime;
            if (w->remaining == 1) { ++w->recover_index; w->unit_mode = 0; }
        }
        w->operation = CO_RECOVER; return 0;
    }
    if (w->operation == CO_REPLAY_CHECK) {
        mpz_gcd(w->divisor, w->temporary[1], c->n);
        if (proper(c, w->divisor)) { curve_finish(c, w->divisor); return 0; }
        if (mpz_cmp(w->divisor, c->n) == 0) { curve_finish(c, NULL); return 0; }
        pair_set(w->replay, w->temporary);
        const b4_curve_unit *u = &w->powers[i->replay_index];
        i->replay_power *= u->prime;
        if (i->replay_power == u->power) { ++i->replay_index; i->replay_power = 1; }
        w->operation = CO_IDLE; i->boundary = 1; return 0;
    }
    if (w->operation == CO_GROUP) {
        size_t at = i->index + w->pending;
        if (!w->block->prepared || w->block->centers[at] != i->center) return -1;
        size_t offset = w->block->offsets[at];
        if (!offset || offset > i->distance) return -1;
        mpz_mul(c->scratch[0], w->giant[0], w->bz[offset]);
        /* The first product is live in scratch[0]; subtract the second
         * product directly, preserving the signed unreduced difference. */
        mpz_submul(c->scratch[0], w->bx[offset], w->giant[1]);
        mpz_mod(w->terms[i->terms_count + w->pending], c->scratch[0], c->n);
        modular_product(c, w->temp_product, w->temp_product, w->terms[i->terms_count + w->pending]); ++w->pending;
        if (i->index + w->pending == w->group_end) {
            i->terms_count += w->pending; i->index = w->group_end;
            mpz_set(w->product, w->temp_product); w->pending = 0;
            w->operation = CO_IDLE; i->boundary = 1;
        }
        return 0;
    }
    if (i->phase == CS_SETUP) return setup_action(c);
    if (i->phase == CS_ONE) {
        int prime = i->powers_count >= w->chunk ? 0 : need_prime(w);
        if (prime == 2) return 0;
        if (prime) {
            size_t count = w->chunk - i->powers_count;
            if (count > w->block->count - i->index) count = w->block->count - i->index;
            if (i->powers_count && i->chain != w->block->chain) return -1;
            i->chain = w->block->chain; i->work += count;
            for (size_t j = 0; j < count; ++j) w->powers[i->powers_count++] = w->block->units[i->index++];
            i->boundary = 1; return 0;
        }
        if (!i->powers_count) {
            i->phase = CS_TWO_SETUP; i->left = i->next = i->b1 + 1;
            i->hi = i->b2 + 1; i->index = 0; w->block = NULL; i->boundary = 1; return 0;
        }
        if (i->chain) {
            i->work += 1;
            for (size_t j = 0; j < i->powers_count; ++j) { i->work += w->powers[j].reserve; w->records[j] = w->powers[j].fast; }
            load_pair(c, 0, w->value); mpz_set_ui(c->guard, 1);
            c->records = w->records; c->count = i->powers_count;
            c->record = c->opcode = c->completed = c->done = 0; c->starting = 1; c->mode = 2;
        } else {
            mpz_set_ui(c->scalar, 1);
            for (size_t j = 0; j < i->powers_count; ++j) mpz_mul_ui(c->scalar, c->scalar, w->powers[j].power);
            i->work += mpz_sizeinbase(c->scalar, 2) + 1; math_begin(c, w->value);
        }
        w->operation = CO_FAST; return 0;
    }
    if (i->phase == CS_REPLAY) {
        if (i->replay_index >= i->powers_count) { curve_finish(c, NULL); return 0; }
        uint64_t prime = w->powers[i->replay_index].prime;
        i->work += bits64(prime) + 1; math_ui(c, w->replay, prime); w->operation = CO_REPLAY; return 0;
    }
    if (i->phase == CS_TWO_SETUP) {
        i->center = i->b1 & 1 ? i->b1 : i->b1 - 1;
        i->distance = sqrt64(i->b2);
        if (i->distance > (i->center - 1) / 2) i->distance = (i->center - 1) / 2;
        if (i->distance < 2) { ++i->work; i->phase = CS_TWO; i->boundary = 1; return 0; }
        int status = curve_table(c, (size_t)i->distance + 1);
        if (status) return status;
        i->work += 2; load_pair(c, 0, w->value); doubling(c, 1, 0);
        mpz_set(w->bx[1], c->x[1]); mpz_set(w->bz[1], c->z[1]);
        i->baby_count = 1; i->phase = CS_BABY; i->boundary = 1; return 0;
    }
    if (i->phase == CS_BABY) {
        size_t at = i->baby_count + 1; i->work += 2;
        mpz_set(c->x[0], w->bx[i->baby_count]); mpz_set(c->z[0], w->bz[i->baby_count]);
        if (at == 2) doubling(c, 3, 0);
        else {
            mpz_set(c->x[1], w->bx[1]); mpz_set(c->z[1], w->bz[1]);
            mpz_set(c->x[2], w->bx[at - 2]); mpz_set(c->z[2], w->bz[at - 2]); addition(c, 3, 0, 1, 2);
        }
        mpz_set(w->bx[at], c->x[3]); mpz_set(w->bz[at], c->z[3]);
        i->baby_count = at; if (at == i->distance) i->phase = CS_GIANT;
        i->boundary = 1; return 0;
    }
    if (i->phase == CS_GIANT) {
        i->work += 2 * bits64(i->center); math_ui(c, w->value, i->center);
        w->operation = CO_GIANT_A; return 0;
    }
    if (i->phase == CS_TERMS) {
        if (i->recovery >= i->terms_count) { curve_finish(c, NULL); return 0; }
        ++i->work; mpz_gcd(w->divisor, w->terms[i->recovery], c->n);
        if (proper(c, w->divisor)) curve_finish(c, w->divisor);
        else { ++i->recovery; i->boundary = 1; }
        return 0;
    }
    if (i->phase != CS_TWO) return -1;
    int prime = i->terms_count == w->batch ? 0 : need_prime(w);
    if (prime == 2) return 0;
    if (i->terms_count == w->batch || !prime) {
        if (!i->terms_count) { curve_finish(c, NULL); return 0; }
        unsigned exhausted = i->terms_count != w->batch && !prime;
        ++i->work; mpz_gcd(w->divisor, w->product, c->n);
        if (proper(c, w->divisor)) curve_finish(c, w->divisor);
        else if (mpz_cmp(w->divisor, c->n) == 0) { i->phase = CS_TERMS; i->recovery = 0; }
        else {
            i->terms_count = 0; mpz_set_ui(w->product, 1);
            if (exhausted) curve_finish(c, NULL);
        }
        i->boundary = 1; return 0;
    }
    uint64_t candidate = w->block->units[i->index].prime;
    if (i->distance < 2) {
        i->work += bits64(candidate); math_ui(c, w->value, candidate); w->operation = CO_SMALL; return 0;
    }
    if (candidate > i->center + 2 * i->distance) {
        i->work += 2; load_pair(c, 0, w->giant); load_pair(c, 2, w->previous);
        mpz_set(c->x[1], w->bx[i->distance]); mpz_set(c->z[1], w->bz[i->distance]);
        addition(c, 3, 0, 1, 2); pair_set(w->previous, w->giant); save_pair(c, w->giant, 3);
        i->center += 2 * i->distance; i->boundary = 1; return 0;
    }
    size_t end = i->index;
    while (end < w->block->count && end - i->index < w->batch - i->terms_count &&
           w->block->units[end].prime <= i->center + 2 * i->distance) ++end;
    i->work += 2 * (end - i->index); w->group_end = end; w->pending = 0;
    mpz_set(w->temp_product, w->product); w->operation = CO_GROUP;
    return 0;
}

static unsigned phase_group(unsigned phase) { return phase == CS_SETUP ? 0 : phase <= CS_REPLAY ? 1 : 2; }
static int curve_slice(b4_context *c, uint64_t wall_end, uint64_t cpu_end,
                        unsigned actions, const atomic_uint *cancel) {
    b4_curve_state *w = c->curve;
    if (!w || !w->active) return -1;
    uint64_t cpu_wall = 0;
    unsigned group = phase_group(w->info.phase);
    uint64_t meter = w->phase_clocks ? nanoseconds(NULL, CLOCK_MONOTONIC) : 0;
    int result = 0;
    w->info.boundary = 0;
    for (unsigned action = 0; action < actions && !w->info.done && !w->info.need; ++action) {
        uint64_t performed = w->info.setup_ops + w->info.stage_one_ops + w->info.stage_two_ops;
        if ((cancel && atomic_load(cancel)) || atomic_load(&w->cancelled) ||
            (w->cancel_after_actions && performed >= w->cancel_after_actions)) {
            result = -3; break;
        }
        uint64_t wall = nanoseconds(c, CLOCK_MONOTONIC);
        if (!wall) { result = -1; break; }
        if (wall >= wall_end) break;
        if (!cpu_wall || wall < cpu_wall || wall - cpu_wall >= 1000000u) {
            uint64_t cpu = nanoseconds(c, CLOCK_PROCESS_CPUTIME_ID);
            if (!cpu) { result = -1; break; }
            if (cpu >= cpu_end) break;
            cpu_wall = wall;
        }
        unsigned old = phase_group(w->info.phase);
        if (old == 0) ++w->info.setup_ops;
        else if (old == 1) ++w->info.stage_one_ops;
        else ++w->info.stage_two_ops;
#ifdef B4_COALESCE_UNIT_YIELDS
        unsigned unit_commit;
        unsigned previous_phase = w->info.phase, previous_operation = w->operation;
        size_t previous_terms = w->info.terms_count;
#endif
        result = curve_action(c);
        if (result) break;
#ifdef B4_COALESCE_UNIT_YIELDS
        /* Consume only a successful normal check, never a recovery frontier.
         * All clocks/action limits and the caller's resource polls remain. */
        unit_commit = !w->atomic_actions && !w->info.done && w->info.boundary &&
            ((previous_phase == CS_ONE && previous_operation == CO_FAST_CHECK &&
              w->info.phase == CS_ONE && w->operation == CO_IDLE && !w->info.powers_count) ||
             (previous_phase == CS_TWO && previous_operation == CO_IDLE && previous_terms &&
              w->info.phase == CS_TWO && w->operation == CO_IDLE && !w->info.terms_count)) &&
            mpz_cmp_ui(w->divisor, 1) == 0;
        if (unit_commit) w->info.boundary = 0;
#endif
        if (w->phase_clocks && phase_group(w->info.phase) != group) {
            uint64_t now = nanoseconds(NULL, CLOCK_MONOTONIC);
            w->info.phase_ns[group] += now - meter; meter = now;
            group = phase_group(w->info.phase);
        }
        if (w->info.boundary && (w->atomic_actions ||
            w->info.done || w->info.phase == CS_TERMS ||
            w->info.phase == CS_REPLAY ||
            (w->info.phase == CS_ONE && !w->info.powers_count && w->operation == CO_IDLE) ||
            (old == 2 && w->info.phase == CS_TWO && !w->info.terms_count))) break;
    }
    if (w->phase_clocks) w->info.phase_ns[group] += nanoseconds(NULL, CLOCK_MONOTONIC) - meter;
    if (result) return result;
    return w->info.done ? 1 : w->info.need ? 2 : 0;
}

int b4_curve_info_get(const b4_context *c, b4_curve_info *out) {
    if (!curve_active(c) || !out) return -1;
    *out = c->curve->info;
    out->strict_state = c->curve->operation == CO_STRICT ? c->curve->strict_phase : 0;
    out->operation = c->curve->operation;
    return 0;
}

int b4_curve_power(const b4_context *c, size_t at, unsigned long long *pair) {
    if (!curve_active(c) || !pair || at >= c->curve->info.powers_count) return -1;
    pair[0] = c->curve->powers[at].prime; pair[1] = c->curve->powers[at].power;
    return 0;
}

int b4_curve_get(const b4_context *c, unsigned field, size_t at, b4_atom *out) {
    if (!curve_active(c) || !out || field > 15) return -1;
    b4_curve_state *w = c->curve;
    if (field >= 14 && (!w->info.done || w->info.phase != CS_REPLAY ||
                         w->operation != CO_REPLAY_CHECK)) return -1;
    mpz_srcptr value = NULL;
    if (field == 0) value = c->a24;
    else if (field <= 2) value = w->value[field - 1];
    else if (field <= 4) value = w->replay[field - 3];
    else if (field == 5) value = w->product;
    else if (field == 6) value = w->factor;
    else if (field <= 8) {
        if (!at || at > w->info.baby_count) return -1;
        value = field == 7 ? w->bx[at] : w->bz[at];
    } else if (field == 9) {
        if (at >= w->info.terms_count) return -1;
        value = w->terms[at];
    } else if (field <= 11) value = w->giant[field - 10];
    else if (field <= 13) value = w->previous[field - 12];
    else value = w->temporary[field - 14];
    export_value(out->data, &out->size, value);
    return 0;
}

static int curve_restore(b4_context *c, const b4_curve_info *s,
                      const b4_value *fields, const b4_value *baby,
                      const b4_value *terms, b4_prime_block *powers,
                      b4_prime_block *block) {
    if (!curve_active(c) || !s || !fields || s->done || s->phase > CS_TERMS ||
        s->b1 != c->curve->info.b1 || s->b2 != c->curve->info.b2 ||
        s->powers_count > c->curve->chunk || s->terms_count > c->curve->batch ||
        s->baby_count > CURVE_TABLE_LIMIT || s->replay_index > s->powers_count ||
        s->recovery > s->terms_count || (s->powers_count && !powers) ||
        (s->baby_count && !baby) || (s->terms_count && !terms) ||
        (s->powers_count && (powers->count != s->powers_count || powers->bound != s->b1 || powers->chain != s->chain))) return -1;
    b4_curve_state *w = c->curve;
    mpz_ptr targets[11] = {c->a24, w->value[0], w->value[1], w->replay[0], w->replay[1],
        w->product, w->factor, w->giant[0], w->giant[1], w->previous[0], w->previous[1]};
    for (unsigned j = 0; j < 11; ++j) {
        if (!canonical(fields[j], 128)) return -1;
        import_value(targets[j], fields[j]);
        if (mpz_cmp(targets[j], c->n) >= 0) return -1;
    }
    if (s->phase != CS_SETUP && !mpz_sgn(w->value[0]) && !mpz_sgn(w->value[1])) return -1;
    uint64_t distance = sqrt64(s->b2), center = s->b1 & 1 ? s->b1 : s->b1 - 1;
    if (distance > (center - 1) / 2) distance = (center - 1) / 2;
    if ((s->phase >= CS_BABY && s->distance != distance) || s->baby_count > distance) return -1;
    if (s->baby_count) {
        int status = curve_table(c, (size_t)distance + 1); if (status) return status;
        for (size_t j = 1; j <= s->baby_count; ++j) {
            if (!canonical(baby[2*(j-1)],128) || !canonical(baby[2*(j-1)+1],128)) return -1;
            import_value(w->bx[j],baby[2*(j-1)]); import_value(w->bz[j],baby[2*(j-1)+1]);
            if (mpz_cmp(w->bx[j],c->n)>=0 || mpz_cmp(w->bz[j],c->n)>=0) return -1;
        }
    }
    for (size_t j = 0; j < s->terms_count; ++j) {
        if (!canonical(terms[j],128)) return -1;
        import_value(w->terms[j],terms[j]); if (mpz_cmp(w->terms[j],c->n)>=0) return -1;
    }
    w->info = *s; w->info.need = 0; w->info.work = 0;
    memset(w->info.phase_ns,0,sizeof(w->info.phase_ns));
    w->block = block;
    if (block) {
        if (s->left != block->lo || s->next != block->hi || s->index > block->count) return -1;
        if (w->block_count >= CURVE_BLOCK_LIMIT) return -2;
        w->blocks[w->block_count++] = block; ++block->refs;
    } else if (s->index) return -1;
    if (s->powers_count) {
        for (size_t j = 0; j < s->powers_count; ++j) w->powers[j] = powers->units[j];
        if (powers != block) {
            if (w->block_count >= CURVE_BLOCK_LIMIT) return -2;
            w->blocks[w->block_count++] = powers; ++powers->refs;
        }
    }
    w->operation = CO_IDLE;
    return 0;
}

int b4_curve_restore(b4_context *c, const b4_curve_info *s,
                      const b4_value *fields, const b4_value *baby,
                      const b4_value *terms, b4_prime_block *powers,
                      b4_prime_block *block) {
    int result = curve_restore(c, s, fields, baby, terms, powers, block);
    if (result) b4_abort(c);
    return result;
}
