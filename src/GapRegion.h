#ifndef REFINE_STAR_GAPREGION_H
#define REFINE_STAR_GAPREGION_H

#include <string>
#include <string_view>
#include <algorithm>
#include <cstdlib>
#include <tuple>
#include "Utils.h"

std::pair<std::vector<std::string>, std::vector<std::string>> slice_alignment(const std::vector<std::string> &sequences, int start, int end) {
    std::vector<std::string> blocks;
    std::vector<std::string> blocks_sequences;
    for (const auto &seq : sequences) {
        std::string curr_segment = std::string(seq).substr(start, end - start + 1);
        blocks.push_back(curr_segment);
        std::string no_gap(curr_segment);
        no_gap.erase(remove(no_gap.begin(), no_gap.end(), '-'), no_gap.end());
        blocks_sequences.push_back(no_gap);
    }
    return make_pair(blocks, blocks_sequences);
}

// Function to find gap regions strictly in a sequence
std::vector<std::pair<int, int>> find_gap_regions_strictly(const std::string_view &sequence) {
    std::vector<std::pair<int, int>> regions;
    int start = -1;

    for (size_t i = 0; i < sequence.length(); ++i) {
        if (sequence[i] == '-') {
            if (start == -1) {
                start = i;
            }
        } else {
            if (start != -1) {
                regions.emplace_back(start, i - 1);
                start = -1;
            }
        }
    }

    if (start != -1) {
        regions.emplace_back(start, sequence.length() - 1);
    }

    return regions;
}


// Function to find gap regions roughly in a sequence
std::vector<std::pair<int, int>> find_gap_regions_roughly(const std::string &sequence, int max_non_gap_bases = 1, int min_region_length = 5) {
    std::vector<std::pair<int, int>> regions;
    int start = -1;
    int non_gap_count = 0;

    for (size_t i = 0; i < sequence.length(); ++i) {
        if (sequence[i] == '-') {
            if (start == -1) {
                start = i;
            }
            non_gap_count = 0; // Reset non-gap counter within a gap region
        } else {
            if (start != -1) {
                non_gap_count++;
                if (non_gap_count > max_non_gap_bases) {
                    // Check if the region length is at least min_region_length
                    if ((i - non_gap_count - start) >= min_region_length) {
                        regions.emplace_back(start, i - non_gap_count);
                    }
                    start = -1;
                    non_gap_count = 0;
                }
            }
        }
    }

    // Check if there's an unclosed gap at the end of sequence and its length
    if (start != -1 && (sequence.length() - start) >= min_region_length) {
        regions.emplace_back(start, sequence.length() - 1);
    }

    return regions;
}

// Function to realign block
std::vector<std::string> realign_block(std::string msa, const std::vector<std::string> &ids, const std::vector<std::string> &sequences, int start, int end) {
    std::vector<std::string> block_sequence;

    if (end - start >= 4) {
        auto before_realign_sequence = slice_alignment(sequences, start, end);
//        auto before_realign_sequence_preprocessed = preprocess(before_realign_sequence.first);
        long long sp_before_realign = score(before_realign_sequence.first, 0, before_realign_sequence.first[0].size());

        utils::Fasta tmp_block;
        std::ofstream ofs("tmp.fasta");
        if (!ofs) {
            std::cerr << "Error: cannot open file tmp.fasta" << std::endl;
            std::cerr << "Please make sure the file path is correct and has appropriate permissions." << std::endl;
            exit(1);
        }
        tmp_block.identifications = ids;
        tmp_block.sequences = before_realign_sequence.second;
        tmp_block.write_to(ofs);
        ofs.close();

        if (msa == "mafft") {
            system("mafft tmp.fasta > tmp.aligned 2> /dev/null");
        } else {
            system("halign -o tmp.aligned tmp.fasta 2> /dev/null");

        }
        
        auto after_realign_sequence = read_from("tmp.aligned");
        long long sp_after_realign = score(after_realign_sequence.sequences, 0, after_realign_sequence.sequences[0].size());

        std::cout << "****************************" << std::endl;
        std::cout << "Block length: " << before_realign_sequence.first[0].length() << std::endl;
        if (sp_after_realign > sp_before_realign) {
            std::cout << "SP before: " << sp_before_realign << std::endl;
            std::cout << "SP after: " << sp_after_realign << std::endl;
            block_sequence = after_realign_sequence.sequences;
        } else {
            block_sequence = std::vector<std::string>(before_realign_sequence.first.begin(), before_realign_sequence.first.end());
        }
    } else {
        auto sliced = slice_alignment(sequences, start, end).first;
        block_sequence = std::vector<std::string>(sliced.begin(), sliced.end());
    }

    return block_sequence;
}

std::vector<std::string> realign_block_muscle(const std::vector<std::string> &ids, const std::vector<std::string> &sequences, int start, int end) {
    std::vector<std::string> block_sequence;

    if (end - start >= 4) {
        auto before_realign_sequence = slice_alignment(sequences, start, end);
//        auto before_realign_sequence_preprocessed = preprocess(before_realign_sequence.first);
        long long sp_before_realign = score(before_realign_sequence.first, 0, before_realign_sequence.first[0].size());

        utils::Fasta tmp_block;
        std::ofstream ofs("tmp.fasta");
        if (!ofs) {
            std::cerr << "Error: cannot open file tmp.fasta" << std::endl;
            std::cerr << "Please make sure the file path is correct and has appropriate permissions." << std::endl;
            exit(1);
        }
        tmp_block.identifications = ids;
        tmp_block.sequences = before_realign_sequence.second;
        tmp_block.write_to(ofs);
        ofs.close();


        system("muscle -in tmp.fasta -out tmp.aligned 2> /dev/null");
        auto realigned_part_sequences = read_from("tmp.aligned");
        int realigned_length = realigned_part_sequences.sequences[0].size();

        std::vector<std::string> after_realign_sequence;
        after_realign_sequence.reserve(ids.size());
        for (const auto &curr_id : ids) {
            int index = find_string_index(curr_id, realigned_part_sequences.identifications);
            if (index == -1) {
                after_realign_sequence.push_back(std::string(realigned_length, '-'));
            } else {
                after_realign_sequence.push_back(realigned_part_sequences.sequences[index]);
            }
        }

        long long sp_after_realign = score(after_realign_sequence, 0, after_realign_sequence[0].size());

        std::cout << "****************************" << std::endl;
        std::cout << "Block length: " << before_realign_sequence.first[0].length() << std::endl;
        if (sp_after_realign > sp_before_realign) {
            std::cout << "SP before: " << sp_before_realign << std::endl;
            std::cout << "SP after: " << sp_after_realign << std::endl;
            block_sequence = after_realign_sequence;
        } else {
            block_sequence = std::vector<std::string>(before_realign_sequence.first.begin(), before_realign_sequence.first.end());
        }
    } else {
        auto sliced = slice_alignment(sequences, start, end).first;
        block_sequence = std::vector<std::string>(sliced.begin(), sliced.end());
    }

    return block_sequence;
}

void join_blocks(std::vector<std::string> &final_sequence, const std::vector<std::string> &block_seqs) {
    for (size_t i = 0; i < block_seqs.size(); ++i) {
        final_sequence[i] += block_seqs[i];
    }
}
#endif //REFINE_STAR_GAPREGION_H
