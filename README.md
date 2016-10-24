Simulated Annealing for complete graphs
===

## Description
This code implements Simulated Annealing for MAX-CUT problems on {+1,-1}-weighted complete graphs, 
which is used in the benchmark study in the paper:

[[http://science.sciencemag.org/content/early/2016/10/19/science.aah4243]](__A coherent Ising machine for 2000-node optimization problems__),  
Takahiro Inagaki, Yoshitaka Haribara, Koji Igarashi, Tomohiro Sonobe, Shuhei Tamate, Toshimori Honjo, Alireza Marandi, Peter L. McMahon, Takeshi Umeki, Koji Enbutsu, Osamu Tadanaga, Hirokazu Takenouchi, Kazuyuki Aihara, Ken-ichi Kawarabayashi, Kyo Inoue, Shoko Utsunomiya, and Hiroki Takesue,  
_Science_ 2016.

## Usage
Please compile the flies with make command (OpenMP is needed).
Please run with these options:

	./main.out <input graph> <num threads> <sync steps> [<target energy>]

In the paper, we ran

	./main.out ./WK2000_1.rud 1 1000000 -60278

where `WK2000_1.rud` is the complete graph with edge weight {+1,-1} (uniform distribution) used in the benchmark. Here, the `<sync steps>` is set to be an arbitrary large value to disable multithreading.

Input file format is the weighted edge list:

	n m 
	v_i1 v_j1 w_{i1j1}
	v_i2 v_j2 w_{i2j2}
	...
	v_im v_jm w_{imjm}

where _n_ is the number of vertices and _m_ is the number of edges.

Note that, if you solve other graphs, the number of vertices should be specified in the `main.cpp`:

	#define NV (2000)
   
since the spins and adjacency matrix are implemented with binary coding. 

## License
MIT License

Copyright (c) 2016 Shuhei Tamate, Tomohiro Sonobe, and Yoshitaka Haribara

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


## Authors
 - S. Tamate  
  - National Institute of Informatics, Hitotsubashi 2-1-2, Chiyoda-ku, Tokyo 101-8403, Japan

 - T. Sonobe  
  - National Institute of Informatics, Hitotsubashi 2-1-2, Chiyoda-ku, Tokyo 101-8403, Japan  
  - ERATO Kawarabayashi Large Graph Project, Japan Science and Technology Agency, Hitotsubashi 2-1-2, Chiyoda-ku, Tokyo 101-8403, Japan

 - Y. Haribara ([haribara@sat.t.u-tokyo.ac.jp](mailto:haribara@sat.t.u-tokyo.ac.jp))  
  - Department of Mathematical Informatics, The University of Tokyo, Hongo 7-3-1, Bunkyo-ku, Tokyo 113-8656, Japan  
  - Institute of Industrial Science, The University of Tokyo, Komaba 4-6-1, Meguro-ku, Tokyo 153-8505, Japan
