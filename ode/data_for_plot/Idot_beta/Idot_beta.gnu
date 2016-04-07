

# -----------------------------------------------------------------------------------------

set terminal postscript eps enhanced color font 'Times-Roman,24'
set xtics font "Times Roman, 30" ;
set ytics font "Times Roman, 30" ;
set xlabel font ',40';
set ylabel font ',40';
set key samplen 1.5
set sty da lp;
#------------------------STYLE---------------------------------------------------------------------------


set style line 1  lc rgb "red" lt 5  ; set style line 2  lc rgb "green" lt 5;
set style line 3  lc rgb "blue" lt 5; set style line 4  lc rgb "magenta" lt 5;
set style line 5  lc rgb "black" lt 5; set style line 6  lc rgb "dark-orange" lt 5;
set style line 7  lc rgb "red" lt 2; set style line 8  lc rgb "green" lt 2;
set style line 9  lc rgb "blue" lt 2; set style line 10  lc rgb "magenta" lt 2;
set style line 11  lc rgb "black" lt 2; set style line 12  lc rgb "dark-orange" lt 2;
set style line 13  lc rgb "red" lt 3; set style line 14  lc rgb "green" lt 3;
set style line 15  lc rgb "blue" lt 3; set style line 16  lc rgb "magenta" lt 3;
set style line 17  lc rgb "black" lt 3; set style line 18  lc rgb "dark-orange" lt 3;
set style line 19 lc rgb "red" lt 4; set style line 20  lc rgb "green" lt 4;
set style line 21  lc rgb "blue" lt 4; set style line 22  lc rgb "magenta" lt 4;
set style line 23  lc rgb "black" lt 4; set style line 24  lc rgb "dark-orange" lt 4;

# -----------------------------------------------------------------------------------------


set term wxt
set autoscale
set pointsize 1;
#set format "%.1f"

b='2';
r='10';
d='0';d1='1';d2='0.5'
k='5';
a0='0';

set key t r;
set xra[0:300]
file0='Idot_beta_r'.r.'d'.d1.'a'.a0.'k'.k.'_si.dat';
file1='Idot_beta_r'.r.'dba'.a0.'k'.k.'_si.dat';


set xlabel ' -{/Symbol b}'
set ylabel '-dI/I'

#set terminal postscript eps enhanced color font 'Times-Roman,24';out_file='../images/Idot_beta_r'.r.'k'.k.'_si.eps';set output out_file;
set key b r;
set xra[0:20]
plot file0 u (-log($1)):(-log($2)) tit '{/Symbol D}=1' ls 13,\
     file1 u (-log($1)):(-log($2)) tit '{/Symbol D}={/Symbol b} ' ls 15
