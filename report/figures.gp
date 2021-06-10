#
# set terminal epslatex input color solid


# files
ics(z,l)=sprintf("../data/output/collected_ics.z-%+i.l-%i.txt", z, l)
dcs(z,l,en)=sprintf("../data/output/dcs.z-%+i.l-%i.en-%.1f.txt", z, l, en)


# parameter set functions
# X(k) = X_{k}, k = 1, ..., n_{X}

# zproj set (-1, +1)
n_z=2
z(k) = 2*k - 3

# lmax set (linear)
n_l=10
l(k) = k - 1

# energy set (quadratic)
eni=0.1
enf=50.0
n_en=20
en(k) = eni + ((k**2) - 1)*((enf - eni)/((n_en**2) - 1))


# common settings
set key inside right top box
set grid xtics ytics
set palette defined (0 "blue" , 1 "red")
unset colorbox


# figure: ics for each zproj
do for [i=1:n_z] {
   l = 5
   z = z(i)

   set output sprintf("figure_ics_z-%+i_l-%i.tex", z, l)

   set title sprintf("Partial ICS [$z_{\rm{proj}} = %i$]", z)
   set xlabel "Energy [eV]"
   set ylabel "ICS [$a_{0}^{2}$]"

   set xrange [0:50]
   set yrange [0.00001:*]
   set logscale y

   plot for [j=3:3+l] \
        ics(z,l) using 1:j title sprintf("$l = %i$", j-3) \
        with lines palette frac ((j-3)/(1.0*(l+1)))

   pause -1

}
