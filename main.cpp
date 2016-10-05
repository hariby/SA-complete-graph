//
//  Simulated Annealing for {+1,-1}-weighted complete graphs
//  Copyright (c) 2016 Shuhei Tamate, Tomohiro Sonobe, and Yoshitaka Haribara
//

//@@ In order to enhance cache hit rate,
//@@ embedding the number of vertices is necessary for MBitSet
#define NV (2000)

//#define MONITOR
#define OUTPUT
#define LIST

#include <cstring>
#include <string>
#include <iostream>
#include <sstream>
#include <list>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <limits>
#include <cmath>
#include <random>
#include <ctime>
#include <climits>
#include <cassert>

//@@
#include <nmmintrin.h> // SSE4.2
#include <omp.h>
#include <sys/time.h>
#include <fcntl.h> // open
#include <unistd.h>

typedef size_t Vertex;

//@@ the sum of weight can be higher than 32 bits because of n = 10^5
//@@ (and the number of edges can be 10^10 > 2^31)
typedef double WeightType;

//@@ clock() cannot be used because of multi-threading
double gettime(){
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_usec / 1e6 + t.tv_sec;
}

//@@ for input buffer
#define BUF_SIZE (1<<20)
char tbuf[BUF_SIZE];
//@@ efficient I/O
class TextInputReader{
private:
  int fd;
  char* buf;
  int buf_size;
  int cur_size;
  int p;
public:
  TextInputReader(int _fd, int bs);
  ~TextInputReader();
  int getLine(char* dst);
};
TextInputReader::TextInputReader(int _fd, int bs){
  fd = _fd;
  buf_size = bs;
  buf = new char[buf_size];
  cur_size = 0;
  p = 0;
}
TextInputReader::~TextInputReader(){
  delete [] buf;
}
int TextInputReader::getLine(char* dst){
  int ret = 0;
  while(1){
    if(p >= cur_size){
      cur_size = read(fd, buf, buf_size);
      if(cur_size <= 0)
        return -1;
      p = 0;
    }
    if(buf[p] == '\n')
      break;
    dst[ret++] = buf[p++];
  }
  p++;
  dst[ret] = '\0';
  return ret;
}

//@@ parser
inline long long int nextll(char*& p){
  if(*p == '\0')
    return -1;
  long long int ret = 0;
  while(!(*p >= '0' && *p <= '9') && *p != '-')p++;
  if(*p == '\0'){
    std::cerr << "no nextll" << std::endl;
    exit(1);
  }
  bool minus = false;
  if(*p == '-'){
    minus = true;
    p++;
  }
  while(*p >= '0' && *p <= '9')
    ret = ret * 10ll + (*p++) - '0';
  return minus ? -ret : ret;
}

//@@ do not change?
typedef long long int mbsb; // MBitSet base type (64 bits)

//@@ Modified BitSet
//@@ bit state -> weight: 0 -> -1, 1 -> 1
class MBitSet{
public:
  MBitSet(){
    n = NV;
    u = sizeof(mbsb) * 8;
    m = (n + u - 1) / u;
    for(size_t i = 0; i < m; i++)a[i] = 0;
  }
  MBitSet(size_t _n){
    n = _n;
    u = sizeof(mbsb) * 8;
    m = (n + u - 1) / u;
    for(size_t i = 0; i < m; i++)a[i] = 0;
  }
  ~MBitSet(){
  }
  MBitSet(const MBitSet& mbs){
    n = mbs.n;
    u = mbs.u;
    m = mbs.m;
    memcpy(a, mbs.a, sizeof(mbsb) * m);
  }
  MBitSet& operator=(const MBitSet& mbs){
    n = mbs.n;
    u = mbs.u;
    m = mbs.m;
    memcpy(a, mbs.a, sizeof(mbsb) * m);
    return (*this);
  }
  int countBits(mbsb val){
    return _mm_popcnt_u64(val);
  }
  size_t size(){return n;}
  bool get(const size_t i){assert(i < n); return (a[i/u]>>(i%u)) & 1;}
  void set(const size_t i, const int bit){
    assert(i < n && (bit == 0 || bit == 1));
    const size_t p = i / u;
    const size_t q = i % u;
    a[p] |= 1ll<<q;
    if(bit == 0)
      a[p] ^= 1ll<<q;
  }
  void flip(const size_t i){a[i/u] ^= 1ll<<(i%u);}
  int xorCount(MBitSet& mbs){
    assert(n == mbs.size());
    int ret = 0;
    for(size_t i = 0; i < m; i++)
      ret += countBits(a[i] ^ mbs.a[i]);
    return ret;
  }
private:
  size_t n; // number of bits
  size_t u; // sizeof(mbsb) * 8
  size_t m; // number of mbsb
  mbsb a[(NV + sizeof(mbsb) * 8 - 1) / (sizeof(mbsb) * 8)];
    // since allocating by new operation is not chache friendly, declaring the array here is suitble (maybe).
};

//@@ shared data for multi-threading
typedef struct SharedData{
public:
  int num_threads;
  int sync_step;
  std::vector<bool> interrupt_flags;
  std::vector<WeightType> best_energies;
  std::vector<MBitSet> best_spins;
  void init(const int _num_threads, const int _sync_step){
    num_threads = _num_threads;
    sync_step = _sync_step;
    interrupt_flags.resize(num_threads, false);
    best_energies.resize(num_threads);
    best_spins.resize(num_threads);
  }
  void clear(const size_t num_site){
    std::fill(interrupt_flags.begin(), interrupt_flags.end(), false);
    std::fill(best_energies.begin(), best_energies.end(), 1<<30);
    MBitSet spin(num_site);
    std::fill(best_spins.begin(), best_spins.end(), spin);
  }
}ShareData;

double calc_beta(size_t step, size_t total_step, double beta_0=4.0)
{
  double beta = beta_0;
  return beta*log(1 + step/static_cast<double>(total_step));
    /* alternative scheduling function */
//  return beta*step/static_cast<double>(total_step);
//  return beta*sqrt(step/static_cast<double>(total_step));
}

std::string RemoveExt(const std::string &filename)
{
    std::string::size_type pos = filename.find_last_of(".");
    if( pos == std::string::npos ){
        return filename;
    }
    return filename.substr(0, pos);
}

std::string RemoveDirName(const std::string &filename)
{
    std::string::size_type pos = filename.find_last_of("/");
    if( pos == std::string::npos ){
        return filename;
    }
    return filename.substr(pos+1);
}

// @@ initial energy calculation with MBitSet
double IsingEnergy(MBitSet& spin, MBitSet* h, size_t num_site){
  double energy = 0;
  for(size_t i = 0; i < num_site; i++){
    bool s = spin.get(i);
    int one = h[i].xorCount(spin) - (s ? 1 : 0);
    int zero = (num_site) - 1 - one;
    if(s){
      energy -= one - zero;
    }else{
      energy += one - zero;
    }
  }
  energy /= 2.0;
  return energy;
}


std::pair< size_t, std::pair< WeightType, MBitSet> > annealing(//WeightedGraph<int>& graph, 
                                           size_t num_flip, const long long int seed, std::ofstream& energy_out,
MBitSet* h,
const size_t num_site, SharedData& sd, int target=INT_MIN, double beta_0=4.0){
  Vertex min_vertex = 0;//graph.get_min_vertex();
  Vertex max_vertex = num_site - 1;//graph.get_max_vertex();
  std::mt19937 engine(seed);
  std::uniform_int_distribution<> rand_vertex(min_vertex, max_vertex);
  std::uniform_int_distribution<> rand_spin(0, 1);
  MBitSet spin(num_site);
  for (size_t i = 0; i < num_site; ++i){
    spin.set(i, rand_spin(engine));
  }
  WeightType energy = (WeightType)IsingEnergy(spin, h, num_site);
  std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

  WeightType best_energy = energy;
  MBitSet best_spin = spin;
  
  double timer_start, timer_temp;
  timer_start = gettime();//clock();
  timer_temp = timer_start;

  const int thread_id = omp_get_thread_num();
#ifdef MONITOR
  #pragma omp critical
  {
    std::cout << "thread_id=" << thread_id << " seed=" << seed << std::endl;
  }
#endif
  size_t flips = 0;
  for (flips = 0; flips < num_flip && !sd.interrupt_flags[thread_id]; ++flips)
  {
    double beta = calc_beta(flips, num_flip, beta_0);

    // choose vertex
    Vertex v = rand_vertex(engine);

    // using MBitSet
    WeightType w_sum3 = 0;
    const int negative = (h[v].xorCount(spin) - (spin.get(v) ? 1 : 0));
    const int positive = (int)num_site - 1 - negative;
    w_sum3 = negative;
    w_sum3 -= positive;
    WeightType del_energy = 2 * (spin.get(v) ? 1 : -1) * w_sum3;

    // accept or not
    if (del_energy > 0)
    {
      double r_a = uniform_dist(engine);
      if (r_a <= exp(-beta*del_energy))
      {
        spin.flip(v);
        energy += del_energy;
      }
    } else {
      spin.flip(v);
      energy += del_energy;

      if (energy < best_energy)
      {
        best_energy = energy;
        best_spin = spin;
      }
    }

    double now = gettime();

#ifdef OUTPUT
    if (now - timer_temp > 0.0001){
      energy_out << (float)(now-timer_start) << ',' << energy << std::endl;
      timer_temp = gettime();//clock();
    }
#endif
    // finish when time>50ms
    if (now - timer_start > 0.05){break;}

    // finish when reached the target value
    if (best_energy <= target){break;}

    // synchronization with other threads
    if(flips > 0 && (flips+1) % sd.sync_step == 0){
      #pragma omp critical
      {
        if(best_energy < sd.best_energies[thread_id]){
          sd.best_energies[thread_id] = best_energy;
          sd.best_spins[thread_id] = best_spin;
        }
        int best_id = thread_id;
        for(int j = 0; j < sd.num_threads; j++){
          if(sd.best_energies[j] < sd.best_energies[best_id])
            best_id = j;
        }
        if(best_id != thread_id){
          energy = sd.best_energies[best_id];
          spin = sd.best_spins[best_id];
        }
      }
    }

  }

#ifdef MONITOR
  if (best_energy > target && !sd.interrupt_flags[thread_id]){std::cout << "insufficient quality " << "thread_id=" << thread_id << std::endl;}
#endif
  return std::pair< size_t, std::pair<WeightType, MBitSet> >(flips, std::pair<WeightType, MBitSet> (best_energy, spin));
}


int main(int argc, char const *argv[])
{
  if (argc < 4)
  {
    std::cout << "Usage: " << argv[0] << " <input graph> <num threads> <sync steps>[<target energy>]\n";
    exit(EXIT_FAILURE);
  }
  std::string filename(argv[1]);

  //@@ setup
  const double setup_start_time = gettime();
  const int num_threads = std::stoi(argv[2]);
  printf("num_threads=%d\n", num_threads);
  omp_set_num_threads(num_threads);
  const int sync_steps = std::stoi(argv[3]);
  printf("sync_steps=%d\n", sync_steps);
  SharedData sd;
  sd.init(num_threads, sync_steps);
  int target = INT_MIN;
  if(argc > 4)
    target = std::stoi(argv[4]);
  printf("target=%d\n", target);
  // read
  int input_fd = open(filename.c_str(), O_RDONLY);
  if(input_fd < 0){
    fprintf(stderr, "cannot open the input file: %s\n", filename.c_str());
    return 1;
  }
  TextInputReader tir(input_fd, 1ll<<27);
  tir.getLine(tbuf);
  char* tbp = tbuf;
  long long int nv = nextll(tbp);
  long long int ne = nextll(tbp);
  printf("nv=%lld ne=%lld\n", nv, ne);
  size_t num_vertices = nv;
  WeightType weight_sum = 0;
  MBitSet* h = new MBitSet[num_vertices];
  
  long long int src, dst;
  int w;
  
  while(tir.getLine(tbuf) >= 0){
    if(tbuf[0] == '#' || tbuf[0] == '%' || tbuf[0] == '*')
      continue;
    tbp = tbuf;
    src = nextll(tbp);
    dst = nextll(tbp);
    w = (int)nextll(tbp);
    src--; dst--;
    assert(src < (int)num_vertices && dst < (int)num_vertices);
    
    const int b = w == 1 ? 1 : 0;
    h[src].set(dst, b);
    h[dst].set(src, b);
    weight_sum += w;
  }
  weight_sum = -weight_sum; //Invert
  
  close(input_fd);
  printf("setup time: %.8f\n", gettime() - setup_start_time);

  int MCS = 200;
  size_t num_flip = MCS * num_vertices;

  double beta_0 = 4.; //.5;

  size_t num_try = 100; //100;
  size_t sum_flip = 0;
  WeightType sum_cut = 0;
  WeightType best_cut = 0;

  //IntVector best_spin(num_vertices);
  MBitSet best_spin(num_vertices);

  float sum_time = 0.;
  float worst_time = 0;
#ifdef LIST
  std::cout << "flip,energy,cut" << std::endl;
#endif
  for (size_t i = 0; i < num_try; ++i){
    sd.clear(num_vertices);
    #pragma omp parallel
    {
      const int thread_id = omp_get_thread_num();
      const long long int seed = i * num_threads + thread_id;
      //clock_t start, end;
      double start, end;
      std::ofstream energy_out;

      std::string energyoutfile("energy_" + RemoveExt(RemoveDirName(filename)) + "_trial" + std::to_string(i) + "_sync-step" + std::to_string(sync_steps) + "_thread" + std::to_string(thread_id) +  ".csv");
#ifdef OUTPUT
      energy_out.open(energyoutfile);
#endif

      start = gettime();//clock();

      std::pair< size_t, std::pair<WeightType, MBitSet> > flip_energy_spin;

      flip_energy_spin = annealing(num_flip, seed, energy_out, h, num_vertices, sd, target, beta_0); // target = -60278 for K_2000 with seed 20001

      //@@ only the first finished thread can enter the following block
      #pragma omp single
      {
#ifdef MONITOR
        printf("thread_id=%d finished\n", thread_id);
#endif
        //@@ set the interruption flags for each thread
        for(int i = 0; i < num_threads; i++)
          sd.interrupt_flags[i] = true;
        size_t flip = flip_energy_spin.first;
        WeightType energy = flip_energy_spin.second.first;
        MBitSet& spin = flip_energy_spin.second.second;
        WeightType cut = -(weight_sum + energy)/2;

#ifdef MONITOR
        printf("flip: %llu\n", (unsigned long long int)flip);
        printf("energy: %.10f\n", (double)energy);
        printf("cut: %.10f\n", (double)cut);
#endif
#ifdef LIST
        std::cout << flip << ',' << (int)energy << ',' << (int)cut << std::endl;
#endif
        if (cut > best_cut){
          best_cut = cut;
          best_spin = spin;
        }
        sum_flip += flip;
        sum_cut += cut;
        end = gettime(); //clock();

        float run_time = (float)(end - start);
        if (run_time > worst_time){
          worst_time = run_time;
        }
        sum_time += run_time;
#ifdef MONITOR
        std::cout << " [Finished in " << run_time << "s]" << std::endl;
#endif
      }
      energy_out.close();
      //@@ wait until all the threads reach here
      #pragma omp barrier
    }
    fflush(stdout);
  }

  std::cout << "\n";
  std::cout << "Results\n";
  std::cout << "-------\n";
  printf("Best: %.10f\n", (double)best_cut);
  printf("Average: %.10f\n", (double)sum_cut/num_try);
  std::cout << "Avg. time: " << (double)sum_time/num_try << std::endl;
  printf("Avg. flip: %.10f\n", (double)sum_flip/num_try);


  std::cout << MCS << ',' << beta_0 << ',' << best_cut << ',' << static_cast<double>(sum_cut)/num_try << ',' << (double)sum_time/num_try << std::endl;

  printf("all finished: %.8f seconds\n", gettime() - setup_start_time);

}
