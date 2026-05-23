// Quick audit: compare unique index keys under murmur3 vs minimap2 hash64 vs raw k-mer.
// Build: clang++ -O2 -std=c++20 -I. -Iinclude hash_collision_audit.cpp sketch.c -o hash_collision_audit
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <string>
#include <unordered_set>
#include <vector>
#include <xseal/xseal_encoder.hpp>
#include <xseal/xseal_minimizer.hpp>
#include <xseal/xseal_kmer_harvester.hpp>

extern "C" {
#include "bseq.h"
#include "mmpriv.h"
}

static inline uint64_t mm2_hash64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask;
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask;
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

static uint64_t kmer_from_ascii(const char *s, int k)
{
	uint64_t mask = (1ULL << (2 * k)) - 1;
	uint64_t f = 0, r = 0;
	for (int i = 0; i < k; ++i) {
		int c = seq_nt4_table[(uint8_t)s[i]];
		if (c >= 4) return UINT64_MAX;
		f = (f << 2 | c) & mask;
		r = (r >> 2) | ((3ULL ^ c) << (2 * (k - 1)));
	}
	return f < r ? f : r;
}

int main(int argc, char **argv)
{
	const char *fn = argc > 1 ? argv[1] : "data/hg38.fa";
	const int k = argc > 2 ? atoi(argv[2]) : 21;
	const int w = argc > 3 ? atoi(argv[3]) : 11;
	const int max_bases = argc > 4 ? atoi(argv[4]) : 50 * 1000 * 1000;

	mm_bseq_file_t *fp = mm_bseq_open(fn);
	if (!fp) {
		fprintf(stderr, "open fail %s\n", fn);
		return 1;
	}
	mm_bseq1_t *s = nullptr;
	int n = 0;
	std::string seq;
	while ((s = mm_bseq_read(fp, 1, 0, &n)) && n > 0) {
		if (s[0].l_seq > 0) {
			seq.append(s[0].seq, s[0].l_seq);
			if ((int)seq.size() >= max_bases) {
				seq.resize(max_bases);
				break;
			}
		}
		free(s[0].seq);
		free(s[0].name);
		free(s);
	}
	mm_bseq_close(fp);
	fprintf(stderr, "Loaded %zu bases from %s\n", seq.size(), fn);

	const uint64_t mask = (1ULL << (2 * k)) - 1;

	// Native minimizers via mm_sketch
	mm128_v mv = {0, 0, 0};
	mm_sketch(0, seq.data(), (int)seq.size(), w, k, 0, 0, &mv);
	const size_t native_hits = mv.n;
	std::unordered_set<uint64_t> native_keys;
	for (uint64_t i = 0; i < mv.n; ++i)
		native_keys.insert(mv.a[i].x >> 8);
	free(mv.a);

	// XSeal path
	xseal::XSealEncoder enc;
	xseal::XSealEncoderState st;
	std::vector<uint8_t> buf(64 * 1024 * 1024);
	enc.init_state(&st, buf.data(), buf.size() - 4096);
	const char *p = seq.data();
	const char *end = p + seq.size();
	enc.encode_chunk(&st, &p, end);
	enc.flush_state(&st);

	xseal::XSealMinimizer scanner(k, w);
	std::vector<uint32_t> pos(8 * 1024 * 1024);
	size_t nh = scanner.scan<11, true>((const __m256i *)buf.data(), st.total_bases_encoded / 128,
	                                   pos.data(), nullptr);
	std::vector<uint64_t> murmur(nh), raw(nh);
	xseal::xseal_harvest_kmers<true>(buf.data(), pos.data(), nh, k, murmur.data());
	xseal::xseal_harvest_kmers<false>(buf.data(), pos.data(), nh, k, raw.data());

	std::unordered_set<uint64_t> u_murmur, u_raw, u_mm2_from_enc, u_mm2_from_ascii;
	std::unordered_set<uint64_t> u_pos;
	uint64_t pos_dup = 0;
	for (size_t i = 0; i < nh; ++i) {
		u_murmur.insert(murmur[i] & mask);
		u_raw.insert(raw[i]);
		u_mm2_from_enc.insert(mm2_hash64(raw[i], mask));
		u_pos.insert(pos[i]);
		if (pos[i] + k <= seq.size()) {
			uint64_t ka = kmer_from_ascii(seq.data() + pos[i], k);
			if (ka != UINT64_MAX)
				u_mm2_from_ascii.insert(mm2_hash64(ka, mask));
		}
	}

	fprintf(stderr, "\n=== Audit (first %d Mb, k=%d w=%d) ===\n", max_bases / 1000000, k, w);
	fprintf(stderr, "native mm_sketch hits: %zu, unique index keys (x>>8): %zu\n", native_hits,
	        native_keys.size());
	fprintf(stderr, "xseal minimizer hits: %zu, unique positions: %zu\n", nh, u_pos.size());
	fprintf(stderr, "xseal unique murmur3 & mask: %zu\n", u_murmur.size());
	fprintf(stderr, "xseal unique raw k-mer (no hash): %zu\n", u_raw.size());
	fprintf(stderr, "xseal unique mm2_hash64(raw_kmer): %zu\n", u_mm2_from_enc.size());
	fprintf(stderr, "xseal unique mm2_hash64(ascii@pos): %zu\n", u_mm2_from_ascii.size());
	fprintf(stderr, "ratio murmur/raw unique: %.3f\n",
	        u_raw.size() ? (double)u_murmur.size() / u_raw.size() : 0.0);
	return 0;
}
