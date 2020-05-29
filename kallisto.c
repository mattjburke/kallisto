/*
 *[KSD] compile with -DKSD to get corrected code.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <stdint.h>

typedef struct{
	char * kmer_seq; // the kmer sequence
	int * t_labels; // list of labels of transcripts that contain the k-mer
	int from_reads;
} kmer_labeled;

typedef struct{
	int eq_class_id; // concatenation of all of the transcript names in the equivalency class
	int count; // the number of fragments(kmers) with the specified equivalency class
	int * eq_class_labels; // list of labels of transcripts that contain the k-mer
} eq_c_count;

uint32_t kmer_hash(char * kmer_seq);
uint32_t eq_hash(char * eq_class);
void store_name_idx(char * name, int t_name_idx);
char** get_kmers(char* seq, int k, int *size);
void create_htable(char * tscripts_file);
void store_read_counts(char * reads_file, int max_reads);
int * get_equiv_class(char* kmer);
void store_eqiv_classes();
uint32_t murmurhash (const char *key, uint32_t len, uint32_t seed);
char *clean(char *str);
void print_eqc_arr();
void print_t_lengths();
void free_hashTable();
void compare(double *probs);


//EM related
void EM(int *lengths, double eps);
void update_parameters(double *alpha, double *probs, int *lengths, double **transcirpt_prob);
void update_alphas(double *alpha, double **transcirpt_prob);
void update_prob(double *alpha, double *probs, int *lengths);
void update_trans_probs(double *probs, double **transcirpt_prob);
double **alloc_matrix(int n, int r);
double likelihood(double *alpha, int *lengths);

kmer_labeled **hashTable;
eq_c_count *eqc_arr;
int tr_name_chars = 13;
int max_hasht_rows = 1000;
int max_per_row = 10;
int k = 50;
int num_t = 125;
int num_eq_classes;
int max_line_length = 7000;
int num_eqc_found;

char ** t_names;
int * t_lengths;


int main()
{

	char * tscripts_file = "transcriptome_gffread.bin";
	create_htable(tscripts_file);

	char * reads_file = "./simdata/cov_100_readlen_150_paired_FALSE/sample_02.fasta";
	store_read_counts(reads_file, INT_MAX);
	store_eqiv_classes();

	// Freeing the hash table because is not needed for the EM
	free_hashTable();

	printf("Eqivalence classes: \n");
	print_eqc_arr();

	EM(t_lengths, 0.01);


	return 0;
}


// returns the index in the hash table that the k-mer sequence and associated transctipt labels are stored at
// use hash function from https://github.com/jwerle/murmurhash.c
// kallisto uses the same murmurhash3 function
uint32_t kmer_hash(char * kmer_seq)
{
	uint32_t seed = 5;
	uint32_t h_val =  murmurhash(kmer_seq, (uint32_t) strlen(kmer_seq), seed) % max_hasht_rows;

	return h_val;
}

// This method returns a list of the possible kmers of a transcripts
// *seq is the transcript, k is the length of the kmer and size if a
// variable used to know the amout of kmers in the transcript
char** get_kmers(char* seq, int k, int *size)
{
	int num_kmers = strlen(seq) - k;
	char ** kmers = malloc(num_kmers * sizeof(char*));

	for (int i = 0; i < num_kmers; i++){
		kmers[i] = malloc(k * sizeof(char));
		strncpy(kmers[i], seq + i, k);
		*size += 1;
	}

	return kmers;
}


// This method returns the amount of characters in 
// the given *file
int fileSize(char *file){
	FILE *fptr;
	char ch;
	fptr=fopen(file,"rb");

	int count = 0;
	num_t = 0;
	while((ch=fgetc(fptr))!=EOF) {
		count++;
		if (ch == '\n'){
			num_t++;
		}
	}
	/* [KSD, BUG] This will fail if there is no newline on the last line. */
	num_t /= 2; // there are two newlines, 1 for the sequence and one for the name, for every transcript in the file
	printf("number of transcripts = %d\n", num_t);
	fclose(fptr);

	return count;	/* [KSD] This is the number of characters in fasta file,
			 * [KSD] which is a VERY ROUGH way to bound the hash size.
			 */
}


//store transcript names in t_names to translate final output. The idex is used as the name id.
void store_name_at_idx(char * name, int t_name_id){
	t_names[t_name_id] = calloc(tr_name_chars, sizeof(char));
	strcpy(t_names[t_name_id], name);
}

// stores the k-mer seqence and the label of the transcript it is read from at index kmer_hash(kmer_seq) in the hash table, if the sequence or T-label are not already present
void create_htable(char * tscripts_file)
//void create_htable()
{
	max_hasht_rows = fileSize(tscripts_file) * 0.5;	/* [KSD] rough */
	t_names = calloc(num_t, sizeof(char*));
	t_lengths = calloc(num_t, sizeof(int));

	printf("Initializing hash table \n");

	hashTable = calloc(max_hasht_rows, sizeof(kmer_labeled*));
	for(int i = 0; i < max_hasht_rows; i++){
		hashTable[i] = calloc(max_per_row, sizeof(kmer_labeled));
		for(int j = 0; j < max_per_row; j++){
			hashTable[i][j].kmer_seq = calloc(k, sizeof(char));
			hashTable[i][j].from_reads = 0;
			hashTable[i][j].t_labels = calloc(num_t, sizeof(int)); //every transcript could contain this seq
			for(int m = 0; m < num_t; m++){
				hashTable[i][j].t_labels[m] = -1; // -1 means null, there is no label
			}
		}
	}
	printf("Hash table initialized\n");

	FILE * tscripts_filep = fopen(tscripts_file, "rb");
	int t_name_id = 0;
	while(!feof(tscripts_filep)){

		char *name	 = malloc(tr_name_chars * sizeof(char));
		char *seq	 = malloc(max_line_length * sizeof(char));

		fgets(name, max_line_length, tscripts_filep);
		fgets(seq, max_line_length, tscripts_filep);

		// cleaning the name for the hash function
		name = clean(name);

		store_name_at_idx(name, t_name_id);
		int l_t = strlen(seq); // needed for EM algo
		t_lengths[t_name_id] = l_t;

		int size = 0;
		char** kmers = get_kmers(seq, k, &size);

		for (int i = 0; i < size; i++){
			char * kmer = malloc(k * sizeof(char));
			strncpy(kmer, kmers[i], k); // just to make code more readable
			uint32_t h_row = kmer_hash(kmer);
			int h_col = 0;

			// find matching sequence or find empty seqence (empty if kmer_seq[0] is null from calloc)
			while (strcmp(hashTable[h_row][h_col].kmer_seq, kmer) != 0 && (hashTable[h_row][h_col].kmer_seq[0]) ){
				h_col++;

				if (h_col >= max_per_row){
					fprintf(stderr, "Not enough columns in hash table: increase max_per_row. \n");
				}
			}

			strcpy(hashTable[h_row][h_col].kmer_seq, kmer); // extra step if found same seq, rewrites seq

			// add transcript id label at first empty spot in array
			int lab_idx = 0;
			while ((hashTable[h_row][h_col].t_labels[lab_idx] != t_name_id) && (hashTable[h_row][h_col].t_labels[lab_idx] != -1 )){
				lab_idx++;
				if (lab_idx >= num_t){ // should never happen
					fprintf(stderr, "More than max number of transcripts have sequence in common. \n");
					lab_idx--;
					break;
				}
			}
			hashTable[h_row][h_col].t_labels[lab_idx] = t_name_id; // extra step if already found, rewrites same int
			//free(kmer);
		}
		t_name_id++;
		//free(name);
		//free(seq);
		//free(kmers);
	}
	printf("Done reading transcripts and storing kmers in hash table!\n\n");
}

void store_read_counts(char * reads_file, int max_reads)
{
	FILE* reads_filep = fopen(reads_file, "rb");

	printf("Reading reads \n");

	int early_stop = 0;
	while(!feof(reads_filep) && early_stop < max_reads){
		early_stop++;

		char *name	 = malloc(max_line_length * sizeof(char));
		char *seq	 = malloc(max_line_length * sizeof(char));

		fgets(name, max_line_length, reads_filep);
		fgets(seq, max_line_length, reads_filep);

		// For the hash table
		name = clean(name);

		// We are using the 1st kmer of the sequence only to identify which transcript the fragment could have come from, as suggested.
		// This is less accurate than the de brujin graph technique kallisto uses, but a necessary shortcut to finish this project in time.
		char* read_1kmer = malloc(k * sizeof(char));
		strncpy(read_1kmer, seq, k);

		// find read_1kmer in hashTable and increment the from_reads count of it
		uint32_t h_row = kmer_hash(read_1kmer);

		int h_col = 0;
		int k_start = 1;
		while (strcmp(hashTable[h_row][h_col].kmer_seq, read_1kmer) != 0) {
			if (!(hashTable[h_row][h_col].kmer_seq[0])){
				char* next_kmer = malloc(k * sizeof(char));	/* [KSD] Why allocate memory instead of directly copying to read_1kmer? */
				strncpy(next_kmer, seq + k_start, k);
				k_start++;
				strcpy(read_1kmer, next_kmer); //begin cycle again with different read_1kmer
				free(next_kmer);
				h_row = kmer_hash(read_1kmer);
				h_col = -1;
			}
			h_col++;
		}

		if (hashTable[h_row][h_col].kmer_seq[0]){
			hashTable[h_row][h_col].from_reads++; // only increment if sequence was found in table
		}
		//free(name);
		//free(seq);
		//free(read_1kmer);
		//printf("hashTable[%d][%d].from_reads = %d\n", h_row, h_col, hashTable[h_row][h_col].from_reads);
	}
	printf("Done storing read counts! \n\n");
}

// returns the list of T-labels associated with the k-mer sequence from the hash table, including -1s
int * get_equiv_class(char* kmer)
{
	uint32_t h_row = kmer_hash(kmer);
	int h_col = 0;
	while (strcmp(hashTable[h_row][h_col].kmer_seq, kmer) != 0) {
		h_col++;
		if (!(hashTable[h_row][h_col].kmer_seq[0])){
			fprintf(stderr, "kmer sequence %s is not stored in hash table.\n", kmer);
			return NULL;
		}
	}
	// allocate memory to copy transcript ids in eqivalency class from hash table
	int * eq_class = malloc(sizeof(int*));
	for (int lab_idx = 0; lab_idx < num_t; lab_idx++){
		eq_class[lab_idx] = hashTable[h_row][h_col].t_labels[lab_idx];
	}
	return eq_class;
}


int find_eqc_arr_idx(int * t_labels){

	num_eqc_found++;
	int idx = num_eqc_found; // new class, return if t_labels doesn't match any previously found classes

	for (int i = 0; i < num_eqc_found; i++){
		for (int j = 0; j < num_t; j++){
			if (eqc_arr[i].eq_class_labels[j] != t_labels[j]){ // they are not the same class, try next in table
				break; // breaks out of both loops?
			}
			else if (j == num_t - 1){ // they are the same class since break statement was never reached
				idx = i;
				num_eqc_found--;
			}
		}
	}
	return idx;
}

// extracting the counts of each eq class is used in the EM algorithm, so this will make getting those faster
void store_eqiv_classes()
{
	printf("Initializing equivalency class array \n");
	num_eq_classes = 400;
	printf("num_eq_classes = %d\n", num_eq_classes);

	if (num_eq_classes < 0){
		fprintf(stderr, "Too many transcrpts, 2^num_t = %d\n", num_eq_classes);
	}

	eqc_arr = calloc(num_eq_classes, sizeof(eq_c_count));
	for(int i = 0; i < num_eq_classes; i++){
		eqc_arr[i].count = 0;
		eqc_arr[i].eq_class_id = -1;
		eqc_arr[i].eq_class_labels = calloc(num_t, sizeof(int));
		for(int j = 0; j < num_t; j++){
			eqc_arr[i].eq_class_labels[j] = -1; // important for comparing equivalency classes
		}
	}

	for(int i = 0; i < max_hasht_rows; i++){
		for(int j = 0; j < max_per_row; j++){
			if (hashTable[i][j].from_reads > 0){
				int idx = -1;
				idx = find_eqc_arr_idx(hashTable[i][j].t_labels); // find the equivalency class of the kmer
				eqc_arr[idx].count += hashTable[i][j].from_reads;

				if (idx == num_eqc_found){
					eqc_arr[idx].eq_class_id = idx; // changes from -1 if new equivalency class
					for(int m = 0; m < num_t; m++){
						eqc_arr[idx].eq_class_labels[m] = hashTable[i][j].t_labels[m];
					}
				}

			}

		}
	}
}

// This method print all the eqivalence classes that were found
void print_eqc_arr(){
	int col_print = 15;

	for (int i = 0 ; i < num_eqc_found; i++){
		printf("eqc_arr[%2d].eq_class_labels = ", i);
		for (int j = 0; j < col_print; j++){
			//if(eqc_arr[i].eq_class_labels[j] != -1)
			printf("%4d ", eqc_arr[i].eq_class_labels[j]);
		}
		printf("count = %6d \n", eqc_arr[i].count);
	}
}

// This method prints the hash table and it's content
void print_hashTable(){
	int col_print = 15;
	for (int i = 0 ; i < max_hasht_rows; i++){
		for (int j = 0; j < max_per_row; j++){
			if (hashTable[i][j].kmer_seq[0]){
				printf("hashTable[%3d][%3d].kmer_seq = %s\n", i, j, hashTable[i][j].kmer_seq);
				printf("hashTable[%3d][%3d].from_reads = %d\n", i, j, hashTable[i][j].from_reads);
				printf("hashTable[%3d][%3d].t_labels = ", i, j); //every transcript could contain this seq
				for(int m = 0; m < col_print; m++){
					printf("%4d ", hashTable[i][j].t_labels[m]); // -1 means null, there is no label
				}
				printf("\n");
			}
		}
	}
}

// Here we are printing the vector of lengths
void print_t_lengths(){
	for (int i = 0; i <= num_t; i++){
		printf("t_lengths[%d] = %d\n", i, t_lengths[i]);
	}
}


// Cleans the given text by changing the next line symbol for a 
// string end symbol
char *clean(char *str){

	for(int i=0; str[i]!='\0'; i++){
		if(str[i]=='\n')
			str[i]='\0';
	}

	return str;

}


// This is an implementation of fast hash table that was used in the paper for kallisto
// from https://github.com/jwerle/murmurhash.c
uint32_t murmurhash (const char *key, uint32_t len, uint32_t seed) {
	uint32_t c1 = 0xcc9e2d51;
	uint32_t c2 = 0x1b873593;
	uint32_t r1 = 15;
	uint32_t r2 = 13;
	uint32_t m = 5;
	uint32_t n = 0xe6546b64;
	uint32_t h = 0;
	uint32_t k = 0;
	uint8_t *d = (uint8_t *) key; // 32 bit extract from `key'
	const uint32_t *chunks = NULL;
	const uint8_t *tail = NULL; // tail - last 8 bytes
	int i = 0;
	int l = len / 4; // chunk length

	h = seed;

	chunks = (const uint32_t *) (d + l * 4); // body
	tail = (const uint8_t *) (d + l * 4); // last 8 byte chunk of `key'

	// for each 4 byte chunk of `key'
	for (i = -l; i != 0; ++i) {
		// next 4 byte chunk of `key'
		k = chunks[i];

		// encode next 4 byte chunk of `key'
		k *= c1;
		k = (k << r1) | (k >> (32 - r1));
		k *= c2;

		// append to hash
		h ^= k;
		h = (h << r2) | (h >> (32 - r2));
		h = h * m + n;
	}

	k = 0;

	// remainder
	switch (len & 3) { // `len % 4'
		case 3: k ^= (tail[2] << 16);
		case 2: k ^= (tail[1] << 8);

		case 1:
			k ^= tail[0];
			k *= c1;
			k = (k << r1) | (k >> (32 - r1));
			k *= c2;
			h ^= k;
	}

	h ^= len;

	h ^= (h >> 16);
	h *= 0x85ebca6b;
	h ^= (h >> 13);
	h *= 0xc2b2ae35;
	h ^= (h >> 16);

	return h;
}


// This is our implementation of the EM algorithm
void EM(int *lengths, double eps){

	// Because the first cell is null in the eqc_arr
	double *probs = calloc(num_t + 1, sizeof(double));
	double *alpha = calloc(num_t + 1,  sizeof(double));

	double llk = 0;
	double prev_llk = 0;
	double **transcirpt_prob = alloc_matrix(num_eqc_found, num_t + 1);

	/* [KSD] Check for failure! */

	// Initialize the probs
	printf("The vaule of num_t is %d\n", num_t + 1);
	for(int i = 0; i <= num_t; i++){
		probs[i] = 1./(num_t + 1);
		alpha[i] = 0;
	}

	double diff;
	do {
		prev_llk = llk;
		update_parameters(alpha, probs, lengths, transcirpt_prob);
		llk = likelihood(alpha, lengths);
		diff = fabs(llk - prev_llk);
fprintf(stderr, "%f %f: %f (%e)\n", alpha[0]/lengths[0], alpha[1]/lengths[1], llk, diff);
	} while(diff > eps);

	// This is the output of the EM
	compare(probs);
	/*
	for(int i = 0; i <= num_t; i++){
		printf("%s \t rho %lf\n", t_names[i], probs[i]);
	}
	*/
	//printf("\n");
}

// This methods update the alphas for each of the transcripts
void update_alphas(double *alpha, double **transcirpt_prob){

	for(int k = 0; k <= num_t; k++){
		double sum = 0;

		for(int i = 1; i < num_eqc_found; i++){
#ifdef KSD	/* [KSD] correction */
			sum += eqc_arr[i].count * transcirpt_prob[i][k];
#else
			sum += transcirpt_prob[i][k];
#endif
		}

		alpha[k] = sum/num_eqc_found;
	}
}

// This updates the probabilities
void update_prob(double *alpha, double *probs, int *lengths){

	double sum = 0;
	for(int j = 0; j <= num_t; j++){
		sum += alpha[j] / lengths[j];
	}

	for(int k = 0; k <= num_t; k++){
		probs[k] = (alpha[k] / lengths[k]) / sum;
	}

}

// Update the transcriptions probabilities
void update_trans_probs(double *probs, double **transcirpt_prob){

	for(int i = 1; i < num_eqc_found; i++){
		int j = 0;
		double sum = 0;

		// Iterate the transcripts
		// Here we are adding the prbability of a specific equivalence class
		while(eqc_arr[i].eq_class_labels[j] != -1){
			sum += probs[eqc_arr[i].eq_class_labels[j]];
			j++;
		}

		for(int k = 0; k <= j; k++)
			transcirpt_prob[i][eqc_arr[i].eq_class_labels[k]] = probs[eqc_arr[i].eq_class_labels[k]] / sum;
	}
}

// Here we are computing the likelihood to determine convergence
double likelihood(double *alpha, int *lengths){
	double llk = 0;

	// Iterate through the equivalence class vector
	for(int e = 1; e < num_eqc_found; e++){
		double sum = 0;
		int t = 0;

		// Iterate through the transcript in a class
		while(eqc_arr[e].eq_class_labels[t] != -1){
			sum += (alpha[eqc_arr[e].eq_class_labels[t]]) / (lengths[eqc_arr[e].eq_class_labels[t]]);
			t++;
		}
		llk += eqc_arr[e].count * log(sum);
	}

	return llk;
}

// From here we call all the update methods
void update_parameters(double *alpha, double *probs, int *lengths, double **transcirpt_prob)
{
	// First, we calculate the probabilities of the trascripts
	update_trans_probs(probs, transcirpt_prob);

	// Second, we update the alphas
	update_alphas(alpha, transcirpt_prob);

	// We update the class probabilities
	update_prob(alpha, probs, lengths);
}

// This method allocate a 2D array of size n and r in the memory
double **alloc_matrix(int r, int c){
	double **mat = calloc(r, sizeof(double*));

	for(int i = 0; i < r; i++)
		mat[i] = calloc(c, sizeof(double));

	return mat;
}

// Here we free all the spaced used by the hash table
void free_hashTable(){
	for(int i = 0; i < max_hasht_rows; i++){
		for(int j = 0; j < max_per_row; j++){
			free(hashTable[i][j].kmer_seq);
			free(hashTable[i][j].t_labels);
		}
		free(hashTable[i]);
	}
	free(hashTable);
}

// This method calculate the mean squared error
void compare(double *probs){


	FILE * tscripts_filep = fopen("real.txt", "r");
	char n[15];
	int d1;
	int *d = malloc(num_t * sizeof(int));
	int tot = 0;

	for(int i = 0; i < num_t; i++){
		fscanf(tscripts_filep, "%s %d %d", n, &d1, &d[i]);
		tot += d[i];
	}

	//Compare the real and the predicted
	double error = 0;
	printf("Labels \t\t Real \t\t Predicted\n");
	for(int i = 0; i < num_t; i++){
		double p = (double) d[i] / tot;
		error += pow((p - probs[i + 1]), 2);
		printf("%s \t %lf \t %lf\n", t_names[i + 1], p, probs[i + 1]);
	}

	printf("The MSE is %lf\n", (error/num_t));


	fclose(tscripts_filep);

}

