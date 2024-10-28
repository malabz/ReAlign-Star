# ReAlign-Star: an optimized realignment method for multiple sequence alignment, targeting star algorithm tools

ReAlign-Star is a tool written in C++17 for realigning the multiple nucleic acid sequence alignment (star algorithm tools). It runs on Linux.

## 🔨Installation and Usage

### 1 Linux/WSL(Windows Subsystem for Linux ) - from the source code

1.Download and Compile the source code. (Make sure your version of gcc >= 9.4.0)
```shell
#1 Download
git clone https://github.com/malabz/ReAlign-Star.git

#2 Open the folder
cd ReAlign-Star

#3 Compile
make

#4 Test ReAlign-Star
./realign_star -h
```

### 2 Usage
```
Usage: ./realign_star -i <input_file> [-o <output_file>] [-w <window_size>] [-l <length>] [-m <msa>]

Options:
  -i <input_file>    (required) Path to the input file containing sequence data.
  -o <output_file>   (optional) Path to the output file for storing results. Default is 'realign_star_result.fasta'.
  -w <window_size>   (optional) Window size for sequence processing. Default is 10.
  -l <length>        (optional) Target length for sequence segments. Default is 5.
  -m <msa>           (optional) MSA tool to use, options are 'halign3', 'mafft', or 'muscle3'. Default is 'mafft'.

Examples:
  ./realign_star -i data.fasta -o results.fasta -w 20 -l 10 -m halign3
  ./realign_star -i data.fasta -m muscle3

Note:
  - The '-i' option is required.
  - The '-m' option only supports 'halign3', 'mafft', and 'muscle3'.
  - If '-w' or '-l' are not provided, default values of 10 and 5 will be used respectively.
```

## 🔬Test dataset and the use case
### 1. Information about the test dataset

Dataset|Sequences Num|Repeats Num|Avg Length|Similarity
:---:|:---:|:---:|:---:|:---:
23S rRNA|500|1|about 3120bp|The average similarity is about 92%
SARS-CoV-2 genome|156|1|about 29000bp|The average similarity is about 99%
SARS-CoV-2 genome|24310|1|about 29000bp|The average similarity is about 99%
16S-like rRNA|100|9|about 1550bp|14 sets of data with different similarities (99%, 98%, 97%, 96%, 95%, 94%, 93%, 92%, 91%, 90%, 85%, 80%, 75%, 70%)
mt-like genome|100|9|about 16000bp|14 sets of data with different similarities (99%, 98%, 97%, 96%, 95%, 94%, 93%, 92%, 91%, 90%, 85%, 80%, 75%, 70%)
SARS-CoV-2-like genome|100|9|about 29000bp|14 sets of data with different similarities (99%, 98%, 97%, 96%, 95%, 94%, 93%, 92%, 91%, 90%, 85%, 80%, 75%, 70%)
CIPRES-128|255|9|about 1550bp|The average similarity is about 80%
CIPRES-256|511|9|about 1550bp|The average similarity is about 80%
CIPRES-512|1023|9|about 1550bp|The average similarity is about 80%
CIPRES-1024|2047|9|about 1550bp|The average similarity is about 80%
CIPRES-2048|4095|9|about 1550bp|The average similarity is about 80%
CIPRES-4096|8191|9|about 1550bp|The average similarity is about 80%

## 📍Reminder
1. Currently ReAlign-Star is **ONLY** available for DNA/RNA. 
3. Please ensure that the sequence ID entered into ReAlign-Star is unique.
4. HAlign3, MAFFT, MUSCLE3 installation is required for the utilization of ReAlign-Star. 

## 🖥️Environment
System|GCC version
:---:|:---:
Linux|GCC 9.4.0
WSL|GCC 9.4.0

## 🔖Citation


## 👋Contacts
The software tools are developed and maintained by 🧑‍🏫[ZOU's lab](http://lab.malab.cn/~zq/en/index.html).

If you find any bug, welcome to contact us on the [issues page](https://github.com/malabz/ReAlign-Star/issues) or email us at 👉[📩](zhai1xiao@gmail.com).

More tools and infomation can visit our [github](https://github.com/malabz).
