#include <iostream>
#include <omp.h>
#include<fstream>
#include <sstream>
#include <vector>
#include <unordered_set>
#include<algorithm>
#include <iomanip>

#include<map>
#include <chrono>  // for high_resolution_clock
using namespace std;

int shared_sid_bid[1000][1000];

std::string revcomp(std::string seq) {
	reverse(seq.begin(), seq.end());
	for (std::size_t i = 0; i < seq.length(); ++i) {
		switch (seq[i]){
		case 'A':
			seq[i] = 'T';
			break;
		case 'C':
			seq[i] = 'G';
			break;
		case 'G':
			seq[i] = 'C';
			break;
		case 'T':
			seq[i] = 'A';
			break;
		}
	}
	return seq;
}

std::map<string, int> generate_bcalm_index(string bcalm_filename, int k, vector<int>& count_bcalm){
	std::map<string, int> kmer_to_bcalm_unitig;
	int uid=0;
	ifstream of_bcalm_another;
	of_bcalm_another.open(bcalm_filename);
	string new_bc_contig;
	 while (std::getline(of_bcalm_another, new_bc_contig)){
			 if(new_bc_contig.at(0)!='>'){
				 count_bcalm.push_back(0);
					for (auto i = 0; i<new_bc_contig.length()-k + 1 ; i++){
							string kmer=new_bc_contig.substr(i, k);
							string rckmer=(revcomp(kmer));
							if(rckmer<kmer)
									kmer = rckmer;
							kmer_to_bcalm_unitig[kmer] = uid;
							count_bcalm[uid]+=1;
					}
					uid++;
			 }
	 }
	 of_bcalm_another.close();

	 cout<<"bcalm index generated"<<endl;
	 return kmer_to_bcalm_unitig;
}


int main(int argc, char **argv){

    if(argc<6){
        cout<<"required args: <k> <cov> <invalid_spades> <invalid_bcalm> <con/uni>"<<endl;
        exit(1);
    }

	for (int i=0; i<1000; i++)
	{
		for (int j=0; j<1000; j++){
			shared_sid_bid[i][j]=0;
		}
	}

    int k=std::stoi(argv[1]);
    int cov=std::stoi(argv[2]);
    //ofstream ofs_gt_spades_bcalm(("gtouch_"+to_string(k)+"_"+to_string(cov)).c_str());

    string con_or_uni = string(argv[5]);
    string output_table_filename = "share_"+con_or_uni+"_"+string(argv[2]);
    ofstream ofs_table(output_table_filename);

		vector<int> count_bcalm;
		vector<int> count_spades;
		std::map<string, int> bcalm_index = generate_bcalm_index(string(argv[4]), k, count_bcalm);


    std::string sp_contig;
		ifstream f_spades(argv[3]);
    int sid=0;
    while (std::getline(f_spades, sp_contig))
    {
        unordered_set<string> spades_single_contig_kmers;
        if(sp_contig.at(0)=='>') continue;

				count_spades.push_back(0);
				for (auto i = 0; i<sp_contig.length()-k + 1 ; i++){
            string kmer=sp_contig.substr(i, k);
            string rckmer=(revcomp(kmer));
            if(rckmer<kmer)
                kmer = rckmer;

						std::map<string,int>::iterator it = bcalm_index.find(kmer);
						if(  it != bcalm_index.end() ){
							int bid = it->second;
							shared_sid_bid[sid][bid]+=1;
						}
            spades_single_contig_kmers.insert(kmer);
        }
      	count_spades[sid]=spades_single_contig_kmers.size();
				sid++;
				std::unordered_set<string>().swap(spades_single_contig_kmers);
    }
    f_spades.close();

		cout<<"spades index done."<<endl;




		cout<<"looping:"<<endl;
		//ofstream

		int bid =count_bcalm.size();

		if(sid>=600000 or bid >=50000){
			cerr<<"array out of bounds"<<endl;
			return EXIT_FAILURE;
		}
		for(int i=0; i<sid; i++){
			for(int j=0; j<bid; j++){
				uint64_t SP_CONTIG_KMER_COUNT =  count_spades[i];
				uint64_t INTERSECT_KMER_COUNT = shared_sid_bid[i][j];
				uint64_t BCALM_KMER_COUNT =  count_bcalm[j];

				ofs_table<<i<<" "<<j<<" ";
				ofs_table<<INTERSECT_KMER_COUNT<<" "<<SP_CONTIG_KMER_COUNT<<" "<<BCALM_KMER_COUNT<<" ";

				double shared=INTERSECT_KMER_COUNT*1.0/SP_CONTIG_KMER_COUNT;
				double shared2=INTERSECT_KMER_COUNT*1.0/BCALM_KMER_COUNT;

				char buffer [50];
				sprintf(buffer, "%0.9f %0.9f \n", shared, shared2);
				ofs_table << string(buffer);
			}
		}

    ofs_table.close();
}
