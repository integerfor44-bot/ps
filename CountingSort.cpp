#include <bits/stdc++.h>
using namespace std;

constexpr int MOD = 998244353;

struct TestCase {
    int n, k;
    vector<int> r;
};

struct Partitions {
    int n;
    vector<int> start, offset, to;
    vector<unsigned char> weight, depth;
    vector<uint64_t> masks;
    unordered_map<uint64_t, int> id;
    int frequency[51]{};

    void generate(int remaining, int minimum, int length,
                  uint64_t mask, int total) {
        if (remaining == 0) {
            int d;

            if (total <= 1) {
                d = 1;
            } else if (length == total) {
                // H(lambda) = [total, 0, ...].
                d = 3;
            } else {
                int parts[51], count = 0;

                for (int x = 1; x <= total; ++x) {
                    if (frequency[x]) {
                        parts[count++] = frequency[x];
                    }
                }

                sort(parts, parts + count);

                uint64_t nextMask = 0;
                for (int i = 0; i < count; ++i) {
                    nextMask |= 1ULL << (parts[i] + i - 1);
                }

                d = 1 + depth[id.at(nextMask)];
            }

            int index = (int)weight.size();
            id.emplace(mask, index);
            masks.push_back(mask);
            weight.push_back((unsigned char)total);
            depth.push_back((unsigned char)d);
            return;
        }

        for (int x = minimum; x <= remaining; ++x) {
            ++frequency[x];

            generate(remaining - x, x, length + 1,
                     mask | (1ULL << (x + length - 1)), total);

            --frequency[x];
        }
    }

    explicit Partitions(int maximum)
        : n(maximum), start(n + 2) {
        int p[51]{};
        p[0] = 1;

        for (int x = 1; x <= n; ++x) {
            for (int s = x; s <= n; ++s) {
                p[s] += p[s - x];
            }
        }

        int states = accumulate(p, p + n + 1, 0);

        weight.reserve(states);
        depth.reserve(states);
        masks.reserve(states);

        id.max_load_factor(0.7f);
        id.reserve(states);

        // State IDs are ordered by partition sum.
        for (int s = 0; s <= n; ++s) {
            start[s] = (int)weight.size();
            generate(s, 1, 0, 0, s);
        }
        start[n + 1] = states;

        offset.resize(states + 1);
        for (int i = 0; i < states; ++i) {
            offset[i + 1] = offset[i] + n - weight[i];
        }

        to.resize(offset.back());

        // Precompute adding a part x to each partition.
        for (int i = 0; i < states; ++i) {
            uint64_t mask = masks[i];
            uint64_t remaining = mask;
            int smaller = 0;

            for (int x = 1; x <= n - weight[i]; ++x) {
                while (remaining &&
                       __builtin_ctzll(remaining) - smaller + 1 < x) {
                    remaining &= remaining - 1;
                    ++smaller;
                }

                // Insert x before all parts >= x.
                int bit = x + smaller - 1;
                uint64_t low = (1ULL << bit) - 1;

                uint64_t nextMask =
                    (mask & low) |
                    ((mask & ~low) << 1) |
                    (1ULL << bit);

                to[offset[i] + x - 1] = id.at(nextMask);
            }
        }

        unordered_map<uint64_t, int>().swap(id);
        vector<uint64_t>().swap(masks);
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int t;
    cin >> t;

    vector<TestCase> tests(t);
    int maximum = 0;

    for (auto &[n, k, r] : tests) {
        cin >> n >> k;
        r.resize(n);

        for (int &x : r) cin >> x;
        maximum = max(maximum, n);
    }

    int choose[51][51]{};

    for (int i = 0; i <= maximum; ++i) {
        choose[i][0] = choose[i][i] = 1;

        for (int j = 1; j < i; ++j) {
            choose[i][j] =
                (choose[i - 1][j - 1] + choose[i - 1][j]) % MOD;
        }
    }

    Partitions partitions(maximum);

    for (const auto &[n, k, r] : tests) {
        int boundCount[51]{};
        for (int x : r) ++boundCount[x];

        int states = partitions.start[n + 1];
        vector<int> dp(states);
        dp[0] = 1;

        int eligible = 0;

        for (int value = n; value >= 1; --value) {
            eligible += boundCount[value];

            // x = 0 is handled by leaving dp[state] unchanged.
            // Positive x increases the state index, so iterate backward.
            for (int state = partitions.start[eligible] - 1;
                 state >= 0; --state) {
                int ways = dp[state];
                if (!ways) continue;

                int available = eligible - partitions.weight[state];
                int edge = partitions.offset[state];

                for (int count = 1; count <= available; ++count) {
                    int destination =
                        partitions.to[edge + count - 1];

                    dp[destination] =
                        (dp[destination] +
                         1LL * ways * choose[available][count]) % MOD;
                }
            }
        }

        vector<int> answer(k + 1);

        auto add = [&](int bucket, int value) {
            if (bucket > k) return;

            answer[bucket] = (answer[bucket] + value) % MOD;
            if (answer[bucket] < 0) answer[bucket] += MOD;
        };

        for (int state = 0; state < states; ++state) {
            add(partitions.depth[state] + 2, dp[state]);
        }

        // These arrays were initially counted in bucket 3:
        // - the zero array;
        // - arrays with exactly one 1 and otherwise zeros.
        int singleOnes = n - boundCount[0];
        int fixedOne = (r[0] >= 1);

        add(3, -1 - singleOnes);
        add(1, 1 + fixedOne);
        add(2, singleOnes - fixedOne);

        for (int p = 1; p <= k; ++p) {
            cout << answer[p] << (p == k ? '\n' : ' ');
        }
    }
}