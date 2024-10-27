#ifndef REFINE_STAR_UTILS_H
#define REFINE_STAR_UTILS_H

#include <fstream>
#include <cmath>
#include <numeric>
#include "Fasta.h"

utils::Fasta read_from(std::string file_path) {
    std::ifstream file(file_path);
    if (!file) {
        std::cerr << "Error: cannot open file " << file_path << std::endl;
        std::cerr << "Please check that the file path is correct and make sure the file exists." << std::endl;
        exit(1);
    }
    utils::Fasta fasta(file);
    file.close();
    return fasta;
}

long long score_column(std::vector<std::string> sequences, unsigned j) {
    static constexpr long long     MISMATCH = -1;
    static constexpr long long        MATCH =  1;
    static constexpr long long GAPEXTENSION = -2;
    static constexpr long long      GAPOPEN = 0;
    unsigned counts[128];
    memset(counts, 0, sizeof(counts));
    for (unsigned i = 0; i != sequences.size(); ++i) {
        char curr_base = sequences[i][j];
        if (curr_base == 'u' || curr_base == 'U') {
            curr_base = 't';
        }
        ++counts[curr_base];
    }

    unsigned const a = counts['a'] + counts['A'];
    unsigned const g = counts['g'] + counts['G'];
    unsigned const c = counts['c'] + counts['C'];
    unsigned const t = counts['t'] + counts['T'];
    unsigned const gap = counts['-'];
    unsigned const n = std::accumulate(counts, counts + 128, 0u) - a - g - c - t - gap;
    return ((a + c) * (g + t) + a * c + g * t) * MISMATCH + ((a * (a - 1) + c * (c - 1) + g * (g - 1) + t * (t - 1)) / 2) * MATCH + ((a + c + g + t + n) * gap) * GAPEXTENSION + ((gap * (gap - 1) / 2) + (n * (n - 1) / 2) + (a + c + g + t) * gap) * GAPOPEN;
}

long long score(std::vector<std::string> sequences, unsigned l, unsigned r) {
    long long s = 0;
    for (unsigned i = l; i != r; ++i)
        s += score_column(sequences, i);
    return s;
}

// Function to remove columns that are all gaps, modifying the input sequences in place
void remove_all_gap_columns(std::vector<std::string>& sequences) {
    if (sequences.empty()) {
        return;  // No operation needed for empty input
    }

    size_t num_sequences = sequences.size();
    size_t sequence_length = sequences[0].size();

    // Create a vector to track which columns are to be kept
    std::vector<bool> keep_column(sequence_length, false);

    // Check each column to see if it contains any non-gap character
    for (size_t j = 0; j < sequence_length; ++j) {
        for (size_t i = 0; i < num_sequences; ++i) {
            if (sequences[i][j] != '-') {
                keep_column[j] = true;
                break;  // No need to check further in this column
            }
        }
    }

    // Filter columns in place
    for (auto& seq : sequences) {
        std::string new_sequence;
        for (size_t j = 0; j < sequence_length; ++j) {
            if (keep_column[j]) {
                new_sequence += seq[j];
            }
        }
        seq = new_sequence;
    }

    // Remove empty sequences if any were left
    sequences.erase(std::remove_if(sequences.begin(), sequences.end(), [](const std::string& seq) {
        return seq.empty();
    }), sequences.end());
}

std::string find_star_sequence(const std::vector<std::string>& sequences) {
    unsigned long long longest_length = 0;
    std::string star_sequence;

    for (const auto &curr : sequences) {
        size_t curr_length = std::count_if(curr.begin(), curr.end(), [](char c) { return c != '-'; });
        if (curr_length > longest_length) {
            star_sequence = curr;
            longest_length = curr_length;
        }
    }

    return star_sequence;
}

int find_string_index(const std::string& a, const std::vector<std::string>& b) {
    // std::find returns an iterator to the found element or b.end() if not found
    auto it = std::find(b.begin(), b.end(), a);

    // Check if the element was found
    if (it != b.end()) {
        // Return the index by subtracting the iterator from b.begin()
        return std::distance(b.begin(), it);
    } else {
        // Return -1 if the element was not found
        return -1;
    }
}

std::vector<int> count_characters_between_dashes(const std::string& input) {
    std::vector<int> result;
    int count = 0;
    bool inside_segment = false;

    for (char c : input) {
        if (c == '-') {
            if (inside_segment && count > 0) {
                result.push_back(count);
            }
            inside_segment = true;
            count = 0;  // Reset the count for the next segment
        } else if (inside_segment) {
            ++count;
        }
    }

    return result;
}

void displayHelp() {
    std::cout << "Usage: program -i <input_file> [-o <output_file>] [-w <window_size>] [-l <length>] [-m <msa>]" << std::endl;
    std::cout << "\nOptions:\n";
    std::cout << "  -i <input_file>    (required) Path to the input file containing sequence data.\n";
    std::cout << "  -o <output_file>   (optional) Path to the output file for storing results. Default is 'realign_star_result.fasta'.\n";
    std::cout << "  -w <window_size>   (optional) Window size for sequence processing. Default is 10.\n";
    std::cout << "  -l <length>        (optional) Target length for sequence segments. Default is 5.\n";
    std::cout << "  -m <msa>           (optional) MSA tool to use, options are 'halign3', 'mafft', or 'muscle3'. Default is 'mafft'.\n";
    std::cout << "\nExamples:\n";
    std::cout << "  ./program -i data.fasta -o results.fasta -w 20 -l 10 -m halign3\n";
    std::cout << "  ./program -i data.fasta -m muscle3\n";
    std::cout << "\nNote:\n";
    std::cout << "  - The '-i' option is required.\n";
    std::cout << "  - The '-m' option only supports 'halign3', 'mafft', and 'muscle3'.\n";
    std::cout << "  - If '-w' or '-l' are not provided, default values of 10 and 5 will be used respectively.\n";
}

#endif //REFINE_STAR_UTILS_H