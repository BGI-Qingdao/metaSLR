#include <iostream>
#include <map>
#include <set>
#include <fstream>
#include <cassert>
#include <ctime>
#include <thread>
#include <stack>
#include <mutex>
#include <vector>
#include <chrono>
#include <getopt.h>
#include "gzstream/gzstream.h"
void logtime() {
    time_t now = time(0);
    char* dt = ctime(&now);
    std::cerr<<dt<<std::endl;
}
//
// load & cache maternal unique kmer & paternal unique kmer
//
std::vector<std::set<std::string>> g_kmers;
int g_K=0;
void load_kmers(const std::string & file,int index){
    assert( g_kmers.size() == index );
    g_kmers.push_back(std::set<std::string>());
    std::ifstream ifs(file);
    std::string line;
    int total_kmer = 0 ;
    if(index==0){
        std::getline(ifs,line);
        g_K = line.size();
        g_kmers[index].insert(line);
        total_kmer++;
    }
    while(!std::getline(ifs,line).eof()){
        g_kmers[index].insert(line);
        total_kmer++;
    }
    std::cerr<<"Recorded "<<total_kmer<<" haplotype "<<index<<" specific "<<g_K<<"-mers\n"; 
}
//
// barcode haplotype relate functions
//
struct BarcodeCache {
    std::map<std::string, std::map<int,int>> barcode_haps;
    void IncrBarcodeHaps(const std::string & barcode , int hap,int incr=1){
        if( barcode_haps[barcode].find(hap) ==  barcode_haps[barcode].end() )
            barcode_haps[barcode][hap] = 0 ;
        barcode_haps[barcode][hap] +=incr ;
    }
    void Add(const BarcodeCache & other){
        for(const auto & pair : other.barcode_haps ) {
            for( const auto & pair1 : pair.second ){
                IncrBarcodeHaps(pair.first,pair1.first,pair1.second);
            }
        }
    }
};

int getHap(const std::string & barcode , const std::map<int,int> & data){
    if( barcode == "0_0_0" || barcode == "0_0" || barcode == "0" )
        return -1;
    int first = 0 , second = 0 , first_index = -1;
    for( int i = 0 ; i< (int)g_kmers.size() ; i++ ){
        if( data.find(i) == data.end() ) 
            continue;
        else if( data.at(i) > first ) {
            second = first ;
            first = data.at(i);
            first_index = i ;
        }
        else if( data.at(i) > second ){
            second = data.at(i);
        }
    }
    if( first > 0 && first_index != -1 && first > second )
        return first_index;
    else 
        return -1 ;
}

int getHapCount(const std::map<int,int> & data , int hap){
    if ( data.find(hap) != data.end() ) return data.at(hap);
    return 0 ;
}

void printBarcodeInfos(const BarcodeCache& g_barcode_haps){
    for(const auto & pair : g_barcode_haps.barcode_haps){
        std::cout<<pair.first;
        const auto & data=pair.second;
        std::cout<<'\t'<<getHap(pair.first,data);
        for( int i = 0 ; i < (int)g_kmers.size() ; i ++ )
            std::cout<<'\t'<<getHapCount(data,i);
        std::cout<<'\t'<<getHapCount(data,-1)<<'\n';
    }
}

//
// kmer relate functions
//
std::map<char,char> g_oppo ;
void InitMap(){
    g_oppo['a']='T';
    g_oppo['A']='T';
    g_oppo['g']='C';
    g_oppo['G']='C';
    g_oppo['c']='G';
    g_oppo['C']='G';
    g_oppo['t']='A';
    g_oppo['T']='A';
}
std::string reverse_complement(const std::string & kmer){
    std::string ret = kmer;
    for(int i = 0 ; i< (int)kmer.size(); i++)
        ret[kmer.size()-i-1] = g_oppo[kmer[i]];
    return ret ;
}

std::string get_cannonical(const std::string & kmer){
   std::string rc = reverse_complement(kmer);
   if (rc < kmer)
      return rc;
   else
      return kmer;
}

//
//reads relate functions
//

// @return : barcode
// A stLFR read's head looks like :
//   @V300017823L1C001R051096800#203_1533_1069/1
//                barcode str :  203_1533_1069
std::string parseName(const std::string & head){
    int s=-1, e=-1;
    for( int i = 0 ; i< (int)head.size(); i++ ){
        if( head[i] == '#' ) s=i;
        if( head[i] == '/' ) e=i;
    }
    return head.substr(s+1,e-s-1);
}

struct MultiThread {
    int t_nums ;
    bool end;
    bool busy;
    void Worker(int index){
        std::pair<std::string,std::string> job;
        while(true){
            locks[index].lock();
            if( caches[index].empty() ){
                busy = false ;
                locks[index].unlock();
                if(end) return ;
                std::this_thread::sleep_for(std::chrono::microseconds(10));
                continue;
            }
            if( ! caches[index].empty() ){
                job=caches[index].top();
                if( caches[index].size() > 10000000 ) busy = true ;
                else if ( caches[index].size() < 300000 ) busy = false ;
                caches[index].pop();
                locks[index].unlock();
                process_reads(job.first,job.second,index);
            } else 
                locks[index].unlock();
        }
    }
    MultiThread(int t_num){
        t_nums = t_num ;
        barcode_caches = new BarcodeCache[t_num];
        locks = new std::mutex[t_num];
        threads = new std::thread*[t_num];
        busy = false;
        end=false;
        for(int i = 0 ; i< t_num ; i++){
            caches.push_back(std::stack< std::pair<std::string,std::string> >());
            threads[i] = new std::thread([this,i](){int index=i ;Worker(index); });
        }
    }
    ~MultiThread(){
        delete [] threads;
        delete [] locks;
        delete [] barcode_caches;
    }
    void process_reads(const std::string & head ,
                         const std::string & seq , int index) {
        std::vector<int> vote;
        for( int i = 0 ; i< (int)g_kmers.size() ; i++ )
            vote.push_back(0);

        for( int i = 0 ; i <(int)seq.size()-g_K+1;i++ ){
            std::string kmer = get_cannonical(seq.substr(i,g_K));
            for( int j = 0 ; j< (int)g_kmers.size() ; j++ ) {
                if( g_kmers[j].find(kmer) != g_kmers[j].end() )
                    vote[j] ++ ;
            }
        }
        int first = 0 , second = 0 , first_index = -1 ;
        for( int i = 0 ; i< (int)g_kmers.size() ; i++ ){
            if( vote[i] > first ) { 
                second = first ;
                first = vote[i];
                first_index = i ;
            }
            else if( vote[i] > second ){
                second = vote[i];
            }
        }
        std::string barcode = parseName(head);
        if (first > 0 && first_index != -1 && first > second ) 
            barcode_caches[index].IncrBarcodeHaps(barcode,first_index);
        else
            barcode_caches[index].IncrBarcodeHaps(barcode,-1);
    }

    void submit(const std::string & head ,const std::string & seq){
        static long index = 0;
        while( busy ) { std::this_thread::sleep_for(std::chrono::seconds(1));}
        index ++ ;
        int id = index % t_nums;
        locks[id].lock();
        caches[id].push(std::make_pair(head,seq));
        locks[id].unlock();
    }
    void wait(){
        for(int i = 0 ; i <t_nums; i++){
            threads[i]->join();
            delete threads[i];
        }
    }
    void collectBarcodes(){
        for(int i = 0 ; i<t_nums ;i++)
            final_data.Add(barcode_caches[i]);
    }
    std::vector<std::stack< std::pair<std::string,std::string> >>  caches;
    std::mutex * locks;
    std::thread ** threads; 
    BarcodeCache * barcode_caches;
    BarcodeCache final_data;
};

void processFastq(const std::string & file,int t_num,BarcodeCache& data){
    std::string head;
    std::string seq;
    std::string tmp;
    MultiThread mt(t_num);
    std::istream *in ;
    bool gz_file = false;
    if( file.size() > 3 ) {
        int end=file.size() ;
        if( file[end-3] == '.' && file[end-2] == 'g' && file[end-1]=='z' ) {
            gz_file = true ;
        }
    }
    if ( gz_file )
        in = new igzstream(file.c_str());
    else 
        in = new std::ifstream(file);
    while(!std::getline(*in,head).eof()){
        std::getline(*in,seq);
        mt.submit(head,seq);
        std::getline(*in,tmp);
        std::getline(*in,tmp);
    }
    mt.end=true;
    mt.wait();
    mt.collectBarcodes();
    data.Add(mt.final_data);
}

void InitAdaptor(){
    std::string r1("CTGTCTCTTATACACATCTTAGGAAGACAAGCACTGACGACATGATCACCAAGGATCGCCATAGTCCATGCTAAAGGACGTCAGGAAGGGCGATCTCAGG");
    std::string r2("TCTGCTGAGTCGAGAACGTCTCTGTGAGCCAAGGAGTTGCTCTGGCGACGGCCACGAAGCTAACAGCCAATCTGCGTAACAGCCAAACCTGAGATCGCCC");
    for(int i = 0 ; i <(int)r1.size()-g_K+1;i++){
        std::string kmer = get_cannonical(r1.substr(i,g_K));
        for ( int j = 0 ; j < (int)g_kmers.size() ; j ++ ){
            if( g_kmers[j].find(kmer) != g_kmers[j].end() ){
                g_kmers[j].erase(kmer);
                std::cerr<<" INFO : erase adaptor kmer from hap "<<j<<" ; kmer="<<kmer<<std::endl;
            }
        }
    }
    for(int i = 0 ; i <(int)r2.size()-g_K+1;i++){
        std::string kmer = get_cannonical(r2.substr(i,g_K));
        for ( int j = 0 ; j < (int)g_kmers.size() ; j ++ ){
            if( g_kmers[j].find(kmer) != g_kmers[j].end() ){
                g_kmers[j].erase(kmer);
                std::cerr<<" INFO : erase adaptor kmer from hap "<<j<<" ; kmer="<<kmer<<std::endl;
            }
        }
    }
}

void printUsage() {
    std::cerr<<"Uasge :\n\tclassify --hap hap0 --hap hap1 [... --hap hapn ] --read read1.fq [--read read2.fq] [--thread t_num]"<<std::endl;
    std::cerr<<"output format: \n\tbarcode haplotype(0/1/2.../n/-1) read_count_hap0 read_count_hap1 ...read_count_hapn read_count_hap-1"<<std::endl;
    std::cerr<<"notice : --read accept file in gzip format , but file must end by \".gz\""<<std::endl;
}
//
// Main function
//

int main(int argc ,char ** argv ){
    static struct option long_options[] = {
        {"hap",  required_argument,  NULL, 'k'},
        {"read", required_argument,  NULL, 'r'},
        {"thread",required_argument, NULL, 't'},
        {"help",  no_argument,       NULL, 'h'},
        {0, 0, 0, 0}
    };
    static char optstring[] = "k:l:r:t:h";
    std::string hap0 , hap1 ;
    std::vector<std::string> haps;
    std::vector<std::string> read;
    int t_num=1;
    while(1){
        int c = getopt_long(argc, argv, optstring, long_options, NULL);
        if (c<0) break;
        switch (c){
            case 'k':
                haps.push_back(std::string(optarg));
                break;
            case 'r':
                read.push_back(std::string(optarg));
                break;
            case 't':
                t_num = atoi(optarg);
                break;
            case 'h':
            default :
                printUsage();
                return -1;
        }
    }
    if( haps.size() < 2 || read.empty() || t_num< 1) {
        printUsage();
        return -1;
    }
    InitMap();
    assert(parseName("VSDSDS#XXX_xxx_s/1")=="XXX_xxx_s");
    assert(get_cannonical("AGCTA")=="AGCTA");
    assert(get_cannonical("TGCTT")=="AAGCA");
    std::cerr<<"__START__"<<std::endl;
    logtime();
    for( int i = 0 ; i < (int)haps.size() ; i++ ) {
        std::cerr<<"__load hap "<<i<<" kmers from file "<<haps[i]<<std::endl;
        load_kmers(haps[i],i);
    }
    InitAdaptor();
    logtime();
    BarcodeCache data;
    for(const auto r : read ){
        std::cerr<<"__process read: "<<r<<std::endl;
        processFastq(r,t_num,data);
        logtime();
        std::cerr<<"__process read done__"<<std::endl;
    }
    std::cerr<<"__print result__"<<std::endl;
    printBarcodeInfos(data);
    logtime();
    std::cerr<<"__END__"<<std::endl;
}