#include <cstdio>
#include <string>
#include <vector>
#include <map>

using namespace std;

map<char, map<int, long long>> counter;

int main(int argc, const char** argv)
{
	string name;
	int min_size = 1;
	long long pos;
	long long homopolymer_start;
	string summary_file = "";
	string output_file = "-";
	vector<string> files;
	bool head_output = false;

	for (int i = 1; i < argc;) {
		if (string(argv[i]) == "-o") {
			output_file = argv[i + 1];
			i += 2;
		} else if (string(argv[i]) == "-s") {
			summary_file = argv[i + 1];
			i += 2;
		} else if (string(argv[i]) == "-n") {
			min_size = atoi(argv[i + 1]);
			i += 2;
		} else {
			files.push_back(argv[i]);
			++i;
		}
	}
	if (files.empty()) {
		fprintf(stderr,
				"\n"
				"Program: Count homopolymer in FASTA files\n"
				"\n"
				"Usage: %s [options] <x.fa>[...]\n"
				"\n"
				"Options:\n"
				"  -n <int>     Minimal length of homopolymer, default: %d\n"
				"  -o <file>    Output the detail, default: %s\n"
				"  -s <file>    Output the summary, default: %s\n"
				"  <x.fa>[...]  Input FASTA files\n"
				"\n", argv[0], min_size, output_file.c_str(), summary_file.c_str());
		return 1;
	}

	FILE *fp_out = NULL;
	if (output_file != "") {
		fp_out = output_file == "-" ? stdout : fopen(output_file.c_str(), "w");
		if (!fp_out) {
			fprintf(stderr, "Can not open file '%s' for writing!\n", output_file.c_str());
			return 1;
		}
	}

	for (size_t i = 0; i < files.size(); ++i) {
		FILE* fp = (files[i] == "-" ? stdin : fopen(files[i].c_str(), "r"));
		if (!fp) {
			fprintf(stderr, "Error: Can not open file '%s'!\n", files[i].c_str());
			return 1;
		}

		char c = fgetc(fp);
		if (c == EOF) {
			fprintf(stderr, "Error: Read file error!\n");
			fclose(fp);
			return 1;
		}
		if (c != '>') {
			fprintf(stderr, "Error: Only support FASTA format!\n");
			fclose(fp);
			return 1;
		}
		name = "";
		while (c > ' ' && !feof(fp)) {
			c = fgetc(fp);
			if (c > ' ')
				name += c;
		}
		while (c != '\n' && !feof(fp)) {
			c = fgetc(fp);
		}
		pos = 0;
		homopolymer_start = 0;

		char last_c = '\0';
		int count = 0;
		if (!head_output) {
			if (fp_out) {
				fprintf(fp_out, "#chr\tpos\tbase\tlength\n");
			}
			head_output = true;
		}
		while (!feof(fp)) {
			c = fgetc(fp);
			if (c >= 'a' && c <= 'z') {
				c -= 'a' - 'A';
			}
			if (c != '>' && c != EOF && c < 'A') {
				continue;
			}
			++pos;

			if (c == last_c) {
				++count;
			} else {
				if (last_c != '\0') {
					if (count >= min_size) {
						if (fp_out) {
							fprintf(fp_out, "%s\t%lld\t%c\t%d\n", name.c_str(), homopolymer_start, last_c, count);
						}
					}
					++counter[last_c][count];
				}
				if (c == '>') {
					name = "";
					while (c > ' ' && !feof(fp)) {
						c = fgetc(fp);
						if (c > ' ')
							name += c;
					}
					while (c != '\n' && !feof(fp)) {
						c = fgetc(fp);
					}
					last_c = '\0';
					pos = 0;
					count = 0;
				} else {
					last_c = c;
					count = 1;
					homopolymer_start = pos;
				}
			}
		}
		fclose(fp);
	}
	if (fp_out && fp_out != stdout) {
		fclose(fp_out);
	}

	if (summary_file != "") {
		fp_out = summary_file == "-" ? stdout : fopen(summary_file.c_str(), "w");
		if (!fp_out) {
			fprintf(stderr, "Error: Can not open file '%s' for writing!\n", summary_file.c_str());
		} else {
			fprintf(fp_out, "#base\tlength\tcount\n");
			for (auto i = counter.begin(); i != counter.end(); ++i) {
				for (auto j = i->second.begin(); j != i->second.end(); ++j) {
					if (j->first >= min_size) {
						fprintf(fp_out, "%c\t%d\t%lld\n", i->first, j->first, j->second);
					}
				}
			}
			if (fp_out != stdout) {
				fclose(fp_out);
			}
		}
	}
	return 0;
}
