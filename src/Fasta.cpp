#include "Fasta.h"

#include <cstring>

utils::Fasta::Fasta(std::istream &is)
{
    _read(is);
}

utils::Fasta::Fasta() {

}

void utils::Fasta::write_to(std::ostream &os, bool with_identification) const
{
    if (with_identification)
        write_to(os, sequences.cbegin(), sequences.cend(), identifications.cbegin());
    else
        write_to(os, sequences.cbegin(), sequences.cend());
}

void utils::Fasta::_read(std::istream &is)
{
    std::string each_line;
    std::string each_sequence;
    size_t longest_length = 0;
    std::string curr = "";
    for (bool flag = false; std::getline(is, each_line); )
    {
        if (each_line.size() == 0)
            continue;

        if (each_line[0] == '>')
        {
            identifications.push_back(each_line.substr(1));
            if (flag) {
                sequences.push_back(std::move(each_sequence));
//                curr = std::move(each_sequence);
//                size_t curr_length = std::count_if(curr.begin(), curr.end(), [](char c) { return c != '-'; });
//                if (curr_length > longest_length) {
//                    longest_length = curr_length;
//                    star_sequence = curr;
//                }
            }
            flag = true;
        }
        else if (flag)
        {
            each_sequence += each_line;
        }
    }

    sequences.push_back(each_sequence);
//    curr = each_sequence;
//    size_t curr_length = count_if(curr.begin(), curr.end(), [](char c) { return c != '-'; });
//    if (curr_length > longest_length) {
//        star_sequence = curr;
//    }
}


void utils::Fasta::cut_and_write(std::ostream &os, const std::string &sequence)
{
    const size_t sequence_length = sequence.size();

    char *cut_sequence = new char[sequence_length + sequence_length / max_line_length + 1];
    size_t des_index = 0;
    for (size_t src_index = 0; src_index < sequence_length; src_index += max_line_length)
    {
        if (src_index) cut_sequence[des_index++] = '\n';

        size_t write_length = sequence_length - src_index;
        if (write_length > max_line_length) write_length = max_line_length;

        memcpy(cut_sequence + des_index, sequence.data() + src_index, write_length);
        des_index += write_length;
    }
    cut_sequence[des_index] = 0;

    os << cut_sequence;
    delete[] cut_sequence;
}

