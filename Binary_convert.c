#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <cmph.h>
#include <stdlib.h>
#include <string.h>


// Prototypes
void write_bin();
void keys_bin();
char **get_kmers(char *seq);


const int STEPSIZE = 100;
//Change this to change the key length
const int kLen = 50;
int max_line_length = 7000;

int main(void){

	write_bin();
	keys_bin();

	return 0;
}

void write_bin(){

	printf("Start writing.\n");

	FILE *pFile = fopen("simdata/transcriptome_gffread.fasta", "r");
	FILE *ptr_myfile = fopen("transcriptome_gffread.bin", "wb");

	/* [KSD] Check for failure to open. */

	int count = 0;
	while(!feof(pFile)){
		char *name	 = malloc(max_line_length * sizeof(char));
		char *seq	 = malloc(max_line_length * sizeof(char));

		fgets(name, max_line_length, pFile);
		fgets(seq, max_line_length, pFile);

		/* [KSD] This is not a binary write.  This file is 
		 * [KSD] identical to the original fasta file:
		 * [KSD] diff simdata/transcriptome_gffread.fasta transcriptome_gffread.bin
		 */
		fprintf(ptr_myfile, name, count);
		fprintf(ptr_myfile, seq, count);

		count++;
	}

	printf("Done writing!\n\n");

	fclose(ptr_myfile);
	fclose(pFile);
}

void keys_bin(){

	printf("Start writing keys.\n");

	FILE *pFile = fopen("transcriptome_gffread.bin" , "rb" );
	FILE *keys = fopen("keys.bin","wb");

	int count = 0;
	while(!feof(pFile)){
		char *name	 = malloc(max_line_length * sizeof(char));
		char *seq	 = malloc(max_line_length * sizeof(char));

		fgets(name, max_line_length, pFile);
		fgets(seq, max_line_length, pFile);

		char **list_k_mers = get_kmers(seq);

		for(int i = 0; i < (strlen(seq) - kLen); i++){

			//printf("%s \n", name);
			//printf("%s \n", list_k_mers[i]);
			fprintf (keys, list_k_mers[i], count++);

			free(list_k_mers[i]);	/* [KSD] added */
		}

		/* [KSD] Serious memory leak!  You need to free
		 * [KSD] list_k_mers after you finish writing them.
		 */
		free(list_k_mers);	/* [KSD] added */

		//printf("%s \n", name);
		//printf("%s \n", seq);
	}

	printf("Done writing keys!\n\n");

	fclose(keys);
	fclose(pFile);

}

char **get_kmers(char *seq){

	int n = strlen(seq);
	char **list = malloc((n - kLen)* sizeof(char*));

	//printf("String %s\n", seq);
	//printf("Substrings:\n");
	for (int i = 0; i < n - kLen; i++)
	{
		//printf("index %d: ", i);
		int j = i + kLen;
		//Print char from current location to ending location
		char *combination = malloc((kLen + 2) * sizeof(char));
		for (int k = i; k < j; k++)
			combination[k - i] = seq[k];

		combination[kLen] = '\n';
		combination[kLen + 1] = '\0';
		list[i] = combination;

	}

	return list;
}
