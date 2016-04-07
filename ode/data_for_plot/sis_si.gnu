

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

b='1';
r='10';
d='2';d1='1.1';d2='0.5'
k='3';
a0='0';

ind(x,A)=(x>A);
ind2(x,A)=(x<A);

file0='dyn_b'.b.'r'.r.'d'.d1.'a'.a0.'k'.k.'_sis.dat';
file2='dyn_b'.b.'r'.r.'d'.d2.'a'.a0.'k'.k.'_sis.dat';
file1='dyn_b'.b.'r'.r.'d'.d.'a'.a0.'k'.k.'_si.dat';

set xlabel ' t'
set ylabel 'log(I)'

#set terminal postscript eps enhanced color font 'Times-Roman,24';out_file='../images/logS_b'.b.'r'.r.'d'.d.'k'.k.'_sis_si.eps';set output out_file;
set key b r;

set k b l;set xra[5:15];set yra[-15:0.1]#-0.00000001]

f(x)=m*x+q
set ylabel 'log(S)'
fit[12:14] f(x) file0 u 1:(log(1.-$3)) via m,q
plot file0 u 1:(log(1.-$3)) tit 'instantanious {/Symbol D}>1' ls 13,\
     file2 u 1:(log(1.-$3)) tit 'instantanious {/Symbol D}<1 ' ls 15,\
     file1 u 1:(log(1.-$3)) tit 'continuous' ls 14,\
     f(x) notit ls 5    

set k b r;
set ylabel 'log(I)'
set yra[-14:0];set xra[0:18]
fit[1:6] f(x) file0 u 1:(log($3)) via m,q
#set terminal postscript eps enhanced color font 'Times-Roman,24';out_file='../images/logI_b'.b.'r'.r.'d'.d.'k'.k.'_sis_si.eps';set output out_file;
plot file0 u 1:(log($3)) tit 'instantanious {/Symbol D}>1' ls 13,\
     file2 u 1:(log($3)) tit 'instantanious {/Symbol D}<1 ' ls 15,\
     file1 u 1:(log($3)) tit 'continuous' ls 14
#     f(x) notit ls 5    

