#include <iostream>
#include <unordered_set>
#include <cstring>
#include <algorithm>
#include <filesystem>
#include <chrono>
#include "Fasta.h"
#include "GapRegion.h"
#include "Utils.h"
#include "Garbage.h"

std::string tmp_folder;

int main(int argc, char **argv) {
    // Check if the program was called with no arguments or with the "-h" option
    if (argc == 1 || (argc == 2 && std::string(argv[1]) == "-h")) {
        displayHelp();
        return 0;
    }
    
    std::string input_file, window, length, msa;
    std::string output_file = "realign_star_result.fasta";


    bool have_window_size = false;
    bool have_length_size = false;
    bool have_msa = false;

    for (int i = 1; i < argc; i += 2) {
        std::string option = argv[i];
        if (i + 1 < argc) {
            std::string value = argv[i + 1];
            if (option == "-i") {
                input_file = value;
            } else if (option == "-o") {
                output_file = value;
            } else if (option == "-w") {
                have_window_size = true;
                window = value;
            } else if (option == "-l") {
                length = value;
                have_length_size = true;
            } else if (option == "-m") {
                msa = value;
                have_msa = true;
            } else {
                std::cerr << "** Unknown option: " << option << std::endl;
                displayHelp();
                return 1;
            }
        } else {
            std::cerr << "** Missing value for option: " << option << std::endl;
            displayHelp();
            return 1;
        }
    }

    if (! have_window_size) {
        window = "10";
    }

    if (! have_length_size) {
        length = "5";
    }

    if (! have_msa) {
        msa = "mafft";
    }

    // Check if required -i option is provided
    if (input_file.empty()) {
        std::cerr << "** Error: -i option is required." << std::endl;
        displayHelp();
        return 1;
    }

    if (msa != "halign3" && msa != "mafft" && msa != "muscle3") {
        std::cerr << "** Error: This MSA tool is not supported." << std::endl;
        displayHelp();
        return 1; 
    }

    // Step 1: Create a random tmp folder using timestamp and random number
    auto timestamp = std::chrono::system_clock::now().time_since_epoch().count();
    int random_number = rand() % 10000;
    tmp_folder = "/tmp/realign_star_" + std::to_string(timestamp) + "_" + std::to_string(random_number);

    // Step 2: Create the temporary folder
    if (!std::filesystem::create_directory(tmp_folder)) {
        std::cerr << "** Error: Failed to create temporary directory." << std::endl;
        return 1;
    }
    std::cout << "Temporary folder created: " << tmp_folder << std::endl;

    // Define the path to profileAlignment.jar in the install directory
    const std::string jar_path = std::string(std::getenv("HOME")) + "/.realign_star/bin/profileAlignment.jar";

    std::string garbage_file = tmp_folder + "/current_bad_sequence.fasta";
    std::string realigned_profile = tmp_folder + "/realigned_profile.fasta";

    std::vector<std::string> final_sequence;
    utils::Fasta alignment = read_from(input_file);

    std::vector<unsigned int> raw_index(alignment.identifications.size());
    std::iota(raw_index.begin(), raw_index.end(), 0);


    //*********** Find garbage sequences - START ***********//
    std::unordered_set<unsigned int> garbage_index = scan_sequences(alignment.sequences, atoi(window.c_str()));
    //*********** Find garbage sequences -  END  ***********//
    if (msa == "muscle3") {
        if (garbage_index.empty()) {
            
            std::string star_sequence = find_star_sequence(alignment.sequences);
            std::cout << "star sequence: " << star_sequence;
            std::cout << std::endl;

            double distance = 0;

            if (alignment.sequences.size() > 1000) {
                distance = 10;
            } else {
                std::vector<int> base_count = count_characters_between_dashes(star_sequence);
                int sum = std::accumulate(base_count.begin(), base_count.end(), 0);
                distance = static_cast<double>(sum) / base_count.size();
                
                if (distance > 10) {
                    distance = 10;
                }
            }

            std::vector<std::pair<int, int>> gap_regions = find_gap_regions_roughly(star_sequence, distance, atoi(length.c_str()));

            std::cout << "Gap regions: ";
            for (const auto &region : gap_regions) {
                std::cout << "(" << region.first << ", " << region.second << ") ";
            }
            std::cout << std::endl;

            if (gap_regions.empty()) {
                std::ofstream ofs(output_file);
                alignment.write_to(ofs);
                ofs.close();
                std::cout << "No bad blocks to realign." << std::endl;
                std::filesystem::remove_all(tmp_folder);
                exit(0);
            }

            if (final_sequence.empty()) {
                final_sequence.resize(alignment.sequences.size(), "");
            }

            if (gap_regions.size() == 1) {
                int curr_start = gap_regions[0].first;
                int curr_end = gap_regions[0].second;
                if (curr_start == 0) {
                    auto realigned_block = realign_block_muscle(alignment.identifications, alignment.sequences, curr_start, curr_end);
                    join_blocks(final_sequence, realigned_block);
                    auto non_realign_region = slice_alignment(alignment.sequences, curr_end + 1, alignment.sequences[0].length() - 1).first;
                    join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                } else {
                    auto non_realign_region = slice_alignment(alignment.sequences, 0, curr_start - 1).first;
                    join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                    auto realigned_block = realign_block_muscle(alignment.identifications, alignment.sequences,curr_start, curr_end);
                    join_blocks(final_sequence, realigned_block);
                    non_realign_region = slice_alignment(alignment.sequences, curr_end + 1, alignment.sequences[0].length() - 1).first;
                    join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                }
                std::ofstream ofs(output_file);
                alignment.sequences = final_sequence;
                alignment.write_to(ofs);
                ofs.close();
                std::filesystem::remove_all(tmp_folder);
                exit(0);
            }

            for (size_t i = 0; i < gap_regions.size(); ++i) {
                int curr_start = gap_regions[i].first;
                int curr_end = gap_regions[i].second;

                if (i == 0) {
                    if (curr_start == 0) {
                        auto realigned_block = realign_block_muscle(alignment.identifications, alignment.sequences, curr_start, curr_end);
                        join_blocks(final_sequence, realigned_block);
                        auto non_realign_region = slice_alignment(alignment.sequences, curr_end + 1, gap_regions[i + 1].first - 1).first;
                        join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                    } else {
                        auto non_realign_region = slice_alignment(alignment.sequences, 0, curr_start - 1).first;
                        join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                        auto realigned_block = realign_block_muscle(alignment.identifications, alignment.sequences,curr_start, curr_end);
                        join_blocks(final_sequence, realigned_block);
                    }
                } else if (i != gap_regions.size() - 1) {
                    if (gap_regions[0].first == 0) {
                        auto realigned_block = realign_block_muscle(alignment.identifications, alignment.sequences,curr_start, curr_end);
                        join_blocks(final_sequence, realigned_block);
                        auto non_realign_region = slice_alignment(alignment.sequences, curr_end + 1, gap_regions[i + 1].first - 1).first;
                        join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                    } else {
                        auto non_realign_region = slice_alignment(alignment.sequences, gap_regions[i - 1].second + 1, curr_start - 1).first;
                        join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                        auto realigned_block = realign_block_muscle(alignment.identifications, alignment.sequences,curr_start, curr_end);
                        join_blocks(final_sequence, realigned_block);
                    }
                } else {
                    if (gap_regions[0].first == 0) {
                        auto realigned_block = realign_block_muscle(alignment.identifications, alignment.sequences,curr_start, curr_end);
                        join_blocks(final_sequence, realigned_block);
                        auto non_realign_region = slice_alignment(alignment.sequences, curr_end + 1, alignment.sequences[0].length() - 1).first;
                        join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                    } else {
                        auto non_realign_region = slice_alignment(alignment.sequences, gap_regions[i - 1].second + 1, curr_start - 1).first;
                        join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                        auto realigned_block = realign_block_muscle(alignment.identifications, alignment.sequences,curr_start, curr_end);
                        join_blocks(final_sequence, realigned_block);
                        non_realign_region = slice_alignment(alignment.sequences, curr_end + 1, alignment.sequences[0].length() - 1).first;
                        join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                    }
                }
            }

            std::ofstream ofs(output_file);
            alignment.sequences = final_sequence;
            alignment.write_to(ofs);
            ofs.close();
            std::filesystem::remove_all(tmp_folder);
            exit(0);
        }
        
        //*********** Find garbage sequences - START ***********//
        std::vector<std::string> garbage_identifications;
        std::vector<std::string> garbage_sequences;
        std::vector<std::string> profile_identifications;
        std::vector<std::string> profile_sequences;

        garbage_identifications.reserve(garbage_index.size());
        garbage_sequences.reserve(garbage_index.size());
        profile_identifications.reserve(alignment.identifications.size() - garbage_index.size());
        profile_sequences.reserve(alignment.identifications.size() - garbage_index.size());

        for (const auto &curr_index : raw_index){
            if (garbage_index.find(curr_index) != garbage_index.end()) {
                garbage_identifications.push_back(alignment.identifications[curr_index]);
                std::string curr_sequence = alignment.sequences[curr_index];
                curr_sequence.erase(std::remove(curr_sequence.begin(), curr_sequence.end(), '-'), curr_sequence.end());
                garbage_sequences.push_back(curr_sequence);
            } else {
                profile_identifications.push_back(alignment.identifications[curr_index]);
                profile_sequences.push_back(alignment.sequences[curr_index]);
            }
        }

        remove_all_gap_columns(profile_sequences);
        //*********** Find garbage sequences - END ***********//

        std::string star_sequence = find_star_sequence(profile_sequences);
        std::cout << "star sequence: " << star_sequence;
        std::cout << std::endl;

        double distance = 0;

        if (alignment.sequences.size() > 1000) {
            distance = 10;
        } else {
            std::vector<int> base_count = count_characters_between_dashes(star_sequence);
            int sum = std::accumulate(base_count.begin(), base_count.end(), 0);
            distance = static_cast<double>(sum) / base_count.size();
                
            if (distance > 10) {
                distance = 10;
            }
        }
        
        std::vector<std::pair<int, int>> gap_regions = find_gap_regions_roughly(star_sequence, distance, atoi(length.c_str()));
        
        std::cout << "Gap regions: ";
        for (const auto &region : gap_regions) {
            std::cout << "(" << region.first << ", " << region.second << ") ";
        }
        std::cout << std::endl;
        
        if (gap_regions.empty()) {
            utils::Fasta profile;
            profile.identifications = profile_identifications;
            profile.sequences = profile_sequences;

            std::ofstream ofs(realigned_profile);
            profile.write_to(ofs);
            ofs.close();
        } else {
            if (gap_regions.size() == 1) {
                int curr_start = gap_regions[0].first;
                int curr_end = gap_regions[0].second;
                if (curr_start == 0) {
                    auto realigned_block = realign_block_muscle(alignment.identifications, alignment.sequences, curr_start, curr_end);
                    join_blocks(final_sequence, realigned_block);
                    auto non_realign_region = slice_alignment(alignment.sequences, curr_end + 1, alignment.sequences[0].length() - 1).first;
                    join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                } else {
                    auto non_realign_region = slice_alignment(alignment.sequences, 0, curr_start - 1).first;
                    join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                    auto realigned_block = realign_block_muscle(alignment.identifications, alignment.sequences,curr_start, curr_end);
                    join_blocks(final_sequence, realigned_block);
                    non_realign_region = slice_alignment(alignment.sequences, curr_end + 1, alignment.sequences[0].length() - 1).first;
                    join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                }
                utils::Fasta profile;
                profile.identifications = profile_identifications;
                profile.sequences = final_sequence;

                std::ofstream ofs(realigned_profile);
                profile.write_to(ofs);
                ofs.close();
            } else {
                if (final_sequence.empty()) {
                    final_sequence.resize(profile_sequences.size(), "");
                }

                for (size_t i = 0; i < gap_regions.size(); ++i) {
                    int curr_start = gap_regions[i].first;
                    int curr_end = gap_regions[i].second;

                    if (i == 0) {
                        if (curr_start == 0) {
                            auto realigned_block = realign_block_muscle(profile_identifications, profile_sequences, curr_start, curr_end);
                            join_blocks(final_sequence, realigned_block);
                            auto non_realign_region = slice_alignment(profile_sequences, curr_end + 1, gap_regions[i + 1].first - 1).first;
                            join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                        } else {
                            auto non_realign_region = slice_alignment(profile_sequences, 0, curr_start - 1).first;
                            join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                            auto realigned_block = realign_block_muscle(profile_identifications, profile_sequences,curr_start, curr_end);
                            join_blocks(final_sequence, realigned_block);
                        }
                    } else if (i != gap_regions.size() - 1) {
                        if (gap_regions[0].first == 0) {
                            auto realigned_block = realign_block_muscle(profile_identifications, profile_sequences,curr_start, curr_end);
                            join_blocks(final_sequence, realigned_block);
                            auto non_realign_region = slice_alignment(profile_sequences, curr_end + 1, gap_regions[i + 1].first - 1).first;
                            join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                        } else {
                            auto non_realign_region = slice_alignment(profile_sequences, gap_regions[i - 1].second + 1, curr_start - 1).first;
                            join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                            auto realigned_block = realign_block_muscle(profile_identifications, profile_sequences,curr_start, curr_end);
                            join_blocks(final_sequence, realigned_block);
                        }
                    } else {
                        if (gap_regions[0].first == 0) {
                            auto realigned_block = realign_block_muscle(profile_identifications, profile_sequences,curr_start, curr_end);
                            join_blocks(final_sequence, realigned_block);
                            auto non_realign_region = slice_alignment(profile_sequences, curr_end + 1, profile_sequences[0].length() - 1).first;
                            join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                        } else {
                            auto non_realign_region = slice_alignment(profile_sequences, gap_regions[i - 1].second + 1, curr_start - 1).first;
                            join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                            auto realigned_block = realign_block_muscle(profile_identifications, profile_sequences,curr_start, curr_end);
                            join_blocks(final_sequence, realigned_block);
                            non_realign_region = slice_alignment(profile_sequences, curr_end + 1, profile_sequences[0].length() - 1).first;
                            join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                        }
                    }
                }

                utils::Fasta profile;
                profile.identifications = profile_identifications;
                profile.sequences = final_sequence;

                std::ofstream ofs(realigned_profile);
                profile.write_to(ofs);
                ofs.close();
            }
        }
        
        utils::Fasta garbages;
        garbages.identifications.resize(1, "");
        garbages.sequences.resize(1, "");

        // Define the path to profileAlignment.jar in the install directory
        const std::string jar_path = std::string(std::getenv("HOME")) + "/.realign_star/bin/profileAlignment.jar";

        // Check if profileAlignment.jar exists
        std::ifstream jar_file(jar_path);
        if (!jar_file.good()) {
            std::cerr << "** Error: profileAlignment.jar not found in " << jar_path << ". Please ensure it is correctly located." << std::endl;
            std::filesystem::remove_all(tmp_folder);
            return 1;
        }
        jar_file.close();

        for(int k = 0; k < garbage_index.size(); k++) {
            garbages.identifications[0] = garbage_identifications[k];
            garbages.sequences[0] = garbage_sequences[k];
            std::ofstream garbage_path(garbage_file);
            garbages.write_to(garbage_path);
            garbage_path.close();

            std::string command_profile_to_seq = "java -jar " + jar_path + " -i " + garbage_file + " " + realigned_profile +" -o " + output_file + " 2> /dev/null";
            system(command_profile_to_seq.c_str());
            std::string command_cp = "cp " + output_file + " " + realigned_profile;
            system(command_cp.c_str());
        }

        std::filesystem::remove_all(tmp_folder);

    } else {
        if (garbage_index.empty()) {
        
        std::string star_sequence = find_star_sequence(alignment.sequences);
        std::cout << "star sequence: " << star_sequence;
        std::cout << std::endl;

        std::vector<int> base_count = count_characters_between_dashes(star_sequence);

        double distance = 0;

        if (alignment.sequences.size() > 1000) {
            distance = 10;
        } else {
            int sum = std::accumulate(base_count.begin(), base_count.end(), 0);
            distance = static_cast<double>(sum) / base_count.size();
            
            if (distance > 10) {
                distance = 10;
            }
        }
        
        std::vector<std::pair<int, int>> gap_regions = find_gap_regions_roughly(star_sequence, distance, atoi(length.c_str()));

        if (gap_regions.empty()) {
            std::ofstream ofs(output_file);
            alignment.write_to(ofs);
            ofs.close();
            std::cout << "No bad blocks to realign." << std::endl;
            std::filesystem::remove_all(tmp_folder);
            exit(0);
        }

        if (final_sequence.empty()) {
            final_sequence.resize(alignment.sequences.size(), "");
        }

        if (gap_regions.size() == 1) {
            int curr_start = gap_regions[0].first;
            int curr_end = gap_regions[0].second;
            if (curr_start == 0) {
                auto realigned_block = realign_block(msa, alignment.identifications, alignment.sequences, curr_start, curr_end);
                join_blocks(final_sequence, realigned_block);
                auto non_realign_region = slice_alignment(alignment.sequences, curr_end + 1, alignment.sequences[0].length() - 1).first;
                join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
            } else {
                auto non_realign_region = slice_alignment(alignment.sequences, 0, curr_start - 1).first;
                join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                auto realigned_block = realign_block(msa, alignment.identifications, alignment.sequences,curr_start, curr_end);
                join_blocks(final_sequence, realigned_block);
                non_realign_region = slice_alignment(alignment.sequences, curr_end + 1, alignment.sequences[0].length() - 1).first;
                join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
            }
            std::ofstream ofs(output_file);
            alignment.sequences = final_sequence;
            alignment.write_to(ofs);
            ofs.close();
            std::filesystem::remove_all(tmp_folder);
            exit(0);
        }

        for (size_t i = 0; i < gap_regions.size(); ++i) {
            int curr_start = gap_regions[i].first;
            int curr_end = gap_regions[i].second;

            if (i == 0) {
                if (curr_start == 0) {
                    auto realigned_block = realign_block(msa, alignment.identifications, alignment.sequences, curr_start, curr_end);
                    join_blocks(final_sequence, realigned_block);
                    auto non_realign_region = slice_alignment(alignment.sequences, curr_end + 1, gap_regions[i + 1].first - 1).first;
                    join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                } else {
                    auto non_realign_region = slice_alignment(alignment.sequences, 0, curr_start - 1).first;
                    join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                    auto realigned_block = realign_block(msa, alignment.identifications, alignment.sequences,curr_start, curr_end);
                    join_blocks(final_sequence, realigned_block);
                }
            } else if (i != gap_regions.size() - 1) {
                if (gap_regions[0].first == 0) {
                    auto realigned_block = realign_block(msa, alignment.identifications, alignment.sequences,curr_start, curr_end);
                    join_blocks(final_sequence, realigned_block);
                    auto non_realign_region = slice_alignment(alignment.sequences, curr_end + 1, gap_regions[i + 1].first - 1).first;
                    join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                } else {
                    auto non_realign_region = slice_alignment(alignment.sequences, gap_regions[i - 1].second + 1, curr_start - 1).first;
                    join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                    auto realigned_block = realign_block(msa, alignment.identifications, alignment.sequences,curr_start, curr_end);
                    join_blocks(final_sequence, realigned_block);
                }
            } else {
                if (gap_regions[0].first == 0) {
                    auto realigned_block = realign_block(msa, alignment.identifications, alignment.sequences,curr_start, curr_end);
                    join_blocks(final_sequence, realigned_block);
                    auto non_realign_region = slice_alignment(alignment.sequences, curr_end + 1, alignment.sequences[0].length() - 1).first;
                    join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                } else {
                    auto non_realign_region = slice_alignment(alignment.sequences, gap_regions[i - 1].second + 1, curr_start - 1).first;
                    join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                    auto realigned_block = realign_block(msa, alignment.identifications, alignment.sequences,curr_start, curr_end);
                    join_blocks(final_sequence, realigned_block);
                    non_realign_region = slice_alignment(alignment.sequences, curr_end + 1, alignment.sequences[0].length() - 1).first;
                    join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                }
            }
        }

        std::ofstream ofs(output_file);
        alignment.sequences = final_sequence;
        alignment.write_to(ofs);
        ofs.close();
        std::filesystem::remove_all(tmp_folder);
        exit(0);
    }
    
    //*********** Find garbage sequences - START ***********//
    std::vector<std::string> garbage_identifications;
    std::vector<std::string> garbage_sequences;
    std::vector<std::string> profile_identifications;
    std::vector<std::string> profile_sequences;

    garbage_identifications.reserve(garbage_index.size());
    garbage_sequences.reserve(garbage_index.size());
    profile_identifications.reserve(alignment.identifications.size() - garbage_index.size());
    profile_sequences.reserve(alignment.identifications.size() - garbage_index.size());

    for (const auto &curr_index : raw_index){
        if (garbage_index.find(curr_index) != garbage_index.end()) {
            garbage_identifications.push_back(alignment.identifications[curr_index]);
            std::string curr_sequence = alignment.sequences[curr_index];
            curr_sequence.erase(std::remove(curr_sequence.begin(), curr_sequence.end(), '-'), curr_sequence.end());
            garbage_sequences.push_back(curr_sequence);
        } else {
            profile_identifications.push_back(alignment.identifications[curr_index]);
            profile_sequences.push_back(alignment.sequences[curr_index]);
        }
    }

    remove_all_gap_columns(profile_sequences);
    //*********** Find garbage sequences - END ***********//

    std::string star_sequence = find_star_sequence(profile_sequences);
    std::cout << "star sequence: " << star_sequence;
    std::cout << std::endl;

    std::vector<int> base_count = count_characters_between_dashes(star_sequence);

    double distance = 0;
    
    if (alignment.sequences.size() > 1000) {
            distance = 10;
    } else {
        int sum = std::accumulate(base_count.begin(), base_count.end(), 0);
        distance = static_cast<double>(sum) / base_count.size();
            
        if (distance > 10) {
            distance = 10;
        }
    }

    std::vector<std::pair<int, int>> gap_regions = find_gap_regions_roughly(star_sequence, distance, atoi(length.c_str()));
    
    std::cout << "Gap regions: ";
    for (const auto &region : gap_regions) {
        std::cout << "(" << region.first << ", " << region.second << ") ";
    }
    std::cout << std::endl;
    
    if (gap_regions.empty()) {
        utils::Fasta profile;
        profile.identifications = profile_identifications;
        profile.sequences = profile_sequences;

        std::ofstream ofs(realigned_profile);
        profile.write_to(ofs);
        ofs.close();
    } else {
        if (gap_regions.size() == 1) {
            int curr_start = gap_regions[0].first;
            int curr_end = gap_regions[0].second;
            if (curr_start == 0) {
                auto realigned_block = realign_block(msa, alignment.identifications, alignment.sequences, curr_start, curr_end);
                join_blocks(final_sequence, realigned_block);
                auto non_realign_region = slice_alignment(alignment.sequences, curr_end + 1, alignment.sequences[0].length() - 1).first;
                join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
            } else {
                auto non_realign_region = slice_alignment(alignment.sequences, 0, curr_start - 1).first;
                join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                auto realigned_block = realign_block(msa, alignment.identifications, alignment.sequences,curr_start, curr_end);
                join_blocks(final_sequence, realigned_block);
                non_realign_region = slice_alignment(alignment.sequences, curr_end + 1, alignment.sequences[0].length() - 1).first;
                join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
            }
            utils::Fasta profile;
            profile.identifications = profile_identifications;
            profile.sequences = final_sequence;

            std::ofstream ofs(realigned_profile);
            profile.write_to(ofs);
            ofs.close();
        } else {
            if (final_sequence.empty()) {
                final_sequence.resize(profile_sequences.size(), "");
            }

            for (size_t i = 0; i < gap_regions.size(); ++i) {
                int curr_start = gap_regions[i].first;
                int curr_end = gap_regions[i].second;

                if (i == 0) {
                    if (curr_start == 0) {
                        auto realigned_block = realign_block(msa, profile_identifications, profile_sequences, curr_start, curr_end);
                        join_blocks(final_sequence, realigned_block);
                        auto non_realign_region = slice_alignment(profile_sequences, curr_end + 1, gap_regions[i + 1].first - 1).first;
                        join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                    } else {
                        auto non_realign_region = slice_alignment(profile_sequences, 0, curr_start - 1).first;
                        join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                        auto realigned_block = realign_block(msa, profile_identifications, profile_sequences,curr_start, curr_end);
                        join_blocks(final_sequence, realigned_block);
                    }
                } else if (i != gap_regions.size() - 1) {
                    if (gap_regions[0].first == 0) {
                        auto realigned_block = realign_block(msa, profile_identifications, profile_sequences,curr_start, curr_end);
                        join_blocks(final_sequence, realigned_block);
                        auto non_realign_region = slice_alignment(profile_sequences, curr_end + 1, gap_regions[i + 1].first - 1).first;
                        join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                    } else {
                        auto non_realign_region = slice_alignment(profile_sequences, gap_regions[i - 1].second + 1, curr_start - 1).first;
                        join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                        auto realigned_block = realign_block(msa, profile_identifications, profile_sequences,curr_start, curr_end);
                        join_blocks(final_sequence, realigned_block);
                    }
                } else {
                    if (gap_regions[0].first == 0) {
                        auto realigned_block = realign_block(msa, profile_identifications, profile_sequences,curr_start, curr_end);
                        join_blocks(final_sequence, realigned_block);
                        auto non_realign_region = slice_alignment(profile_sequences, curr_end + 1, profile_sequences[0].length() - 1).first;
                        join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                    } else {
                        auto non_realign_region = slice_alignment(profile_sequences, gap_regions[i - 1].second + 1, curr_start - 1).first;
                        join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                        auto realigned_block = realign_block(msa, profile_identifications, profile_sequences,curr_start, curr_end);
                        join_blocks(final_sequence, realigned_block);
                        non_realign_region = slice_alignment(profile_sequences, curr_end + 1, profile_sequences[0].length() - 1).first;
                        join_blocks(final_sequence, std::vector<std::string>(non_realign_region.begin(), non_realign_region.end()));
                    }
                }
            }

            utils::Fasta profile;
            profile.identifications = profile_identifications;
            profile.sequences = final_sequence;

            std::ofstream ofs(realigned_profile);
            profile.write_to(ofs);
            ofs.close();
        }
    }
    

    utils::Fasta garbages;
    garbages.identifications.resize(1, "");
    garbages.sequences.resize(1, "");

    // Define the path to profileAlignment.jar in the install directory
    const std::string jar_path = std::string(std::getenv("HOME")) + "/.realign_star/bin/profileAlignment.jar";

    // Check if profileAlignment.jar exists
    std::ifstream jar_file(jar_path);
    if (!jar_file.good()) {
        std::cerr << "** Error: profileAlignment.jar not found in " << jar_path << ". Please ensure it is correctly located." << std::endl;
        std::filesystem::remove_all(tmp_folder);
        return 1;
    }
    jar_file.close();

    for(int k = 0; k < garbage_index.size(); k++) {
        garbages.identifications[0] = garbage_identifications[k];
        garbages.sequences[0] = garbage_sequences[k];
        std::ofstream garbage_path(garbage_file);
        garbages.write_to(garbage_path);
        garbage_path.close();

        std::string command_profile_to_seq = "java -jar " + jar_path + " -i " + garbage_file + " " + realigned_profile +" -o " + output_file + " 2> /dev/null";
        system(command_profile_to_seq.c_str());
        std::string command_cp = "cp " + output_file + " " + realigned_profile;
        system(command_cp.c_str());
    }

    std::filesystem::remove_all(tmp_folder);

    } 

    return 0;
}
