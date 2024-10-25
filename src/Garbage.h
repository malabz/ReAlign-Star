#ifndef REFINE_STAR_GARBAGE_H
#define REFINE_STAR_GARBAGE_H

#include <optional>
#include <vector>
#include <string>
#include <set>
#include <algorithm> // for std::sort
#include <optional>

std::optional<int> is_single_base_sequence(const std::vector<std::string>& region) {
    std::optional<int> base_index = std::nullopt;

    for (int i = 0; i < region.size(); ++i) {
        const std::string& segment = region[i];
        bool has_base = false;

        // 检查这个序列片段是否包含至少一个非 '-' 的字符
        for (char base : segment) {
            if (base != '-') {
                has_base = true;
                break;
            }
        }

        if (has_base) {
            if (base_index.has_value()) {
                // 找到了多个基序列
                return std::nullopt;
            }
            base_index = i;
        }
    }

    return base_index;
}

std::unordered_set<unsigned int> scan_sequences(const std::vector<std::string>& sequences, int window_length) {
    int sequence_length = sequences[0].size();
    std::set<int> seen_sequences;

    // 滑动窗口遍历序列
    for (int start = 0; start <= sequence_length - window_length; ++start) {
        std::vector<std::string> region;

        for (const auto& seq : sequences) {
            region.push_back(seq.substr(start, window_length));
        }

        std::optional<int> base_index = is_single_base_sequence(region);
        if (base_index.has_value() && seen_sequences.find(base_index.value()) == seen_sequences.end()) {
            seen_sequences.insert(base_index.value());
        }
    }

    std::unordered_set<unsigned int> result(seen_sequences.begin(), seen_sequences.end());
    return result;
}





#endif //REFINE_STAR_GARBAGE_H
