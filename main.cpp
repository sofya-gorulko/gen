#include <algorithm>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

struct gen_reading {
    /**
     Struct for saving one reading of gen.
     */
    std::string gen;
    std::string quality;
    int length;

    gen_reading(const std::string &head1, std::string gen, std::string quality) {
        this->gen = std::move(gen);
        this->quality = std::move(quality);
        length = std::stoi(head1.substr(head1.find('=') + 1));
    }

    void trimming(int pos) {
        gen = gen.substr(0, pos);
        quality = quality.substr(0, pos);
        length = gen.size();
    }
};

std::vector<gen_reading> reader(const std::string &file_name) {
    /**
     * Read fastq-file and save results.
     */
    std::ifstream file(file_name);
    if (!file.is_open()) {
        std::cout << "Cannot open file." << std::endl;
        return {};
    }
    std::string head1, head2, gen, quality;
    std::vector<gen_reading> gens;
    while (std::getline(file, head1)) {
        std::getline(file, gen);
        std::getline(file, head2);
        std::getline(file, quality);
        gens.emplace_back(head1, gen, quality);
    }
    file.close();
    return gens;
}

void print_stats(const std::vector<gen_reading> &gens) {
    /**
     * Calculating some statistics: total quantity, minimum length, maximum length, average length.
     */
    int max = 0;
    long long cnt = 0;
    long long sum = 0;
    int min = 1e7;
    for (auto &gen: gens) {
        max = std::max(max, gen.length);
        cnt++;
        sum += gen.length;
        min = std::min(min, gen.length);
    }
    std::cout << "total quantity: " << cnt << std::endl;
    std::cout << "minimum length: " << min << std::endl;
    std::cout << "maximum length: " << max << std::endl;
    std::cout << "average length: " << sum / cnt << std::endl;
}

double gc_stat(const std::vector<gen_reading> &gens) {
    /**
     * Calculate GC-content.
     */
    long long count = 0;
    long long sum = 0;
    for (auto &gen: gens) {
        sum += gen.length;
        count += std::ranges::count(gen.gen, 'C') +
                std::ranges::count(gen.gen, 'G');
    }
    std::cout << "GC-content: " << count * 100. / sum << std::endl;
    return count * 100. / sum;
}

double average_quality_at_position(const std::vector<gen_reading> &gens, int pos) {
    /**
     * Calculate average quality at position {pos}.
     */
    long long cnt = 0;
    long long sum = 0;
    for (auto &gen: gens) {
        sum += gen.quality[pos] - 33;
        cnt++;
    }
    std::cout << "average quality at position " << pos << ": " << sum / cnt << std::endl;
    return sum * 1. / cnt;
}

std::vector<gen_reading> trimming(const std::vector<gen_reading> &gens, int len, int q) {
    /**
     * Trimming with a sliding window {len} and quality {q}.
     */
    std::vector<gen_reading> res;
    int over_quality = len * q;
    for (auto gen: gens) {
        if (gen.length < len)
            continue;
        long long sum = 0;
        int last = gen.length;
        for (int i = 0; i < len; i++) {
            sum += gen.quality[i] - 33;
        }
        if (sum <= over_quality) {
            continue;
        }
        for (int i = len; i < gen.length; i++) {
            sum += gen.quality[i] - gen.quality[i - len];
            if (sum <= over_quality) {
                last = i;
                while (gen.quality[last - 1] - 33 < q && last > 1) {
                    last--;
                }
                if (last != gen.length) {
                    gen.trimming(last);
                }
                break;
            }
        }
        if (last != 0) {
            res.push_back(gen);
        }
    }
    return res;
}

std::vector<gen_reading> trimming(const std::vector<gen_reading> &gens, int len) {
    /**
     * Trimming by length {len}.
     */
    std::vector<gen_reading> res;
    for (auto &gen: gens) {
        if (gen.length >= len)
            res.push_back(gen);
    }
    return res;
}

int main() {
    /**
     * Code for calculate some usefull informwtion about reading of gen.
     */
    std::vector<gen_reading> gens = reader("reads.fastq");
    print_stats(gens);
    gc_stat(gens);
    average_quality_at_position(gens, 9);
    auto trim1 = trimming(gens, 5, 30);
    std::cout << std::endl;
    std::cout << "cut off readings: " << gens.size() - trim1.size() << std::endl;
    print_stats(trim1);
    auto trim2 = trimming(trim1, 60);
    std::cout << "size after trimming: " << trim2.size() << std::endl;
    return 0;
}
