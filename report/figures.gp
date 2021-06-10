#
set terminal epslatex input color solid


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
set palette defined (0 "blue" , 1 "red")
unset colorbox
set grid xtics ytics
set key \
  top right \
  box opaque \
  samplen 1 spacing 0.6 height +0.6


# figure: ics for each zproj
do for [i=1:n_z] {
   z = z(i)
   l = 5

   set output sprintf("figure_ics_z-%+i_l-%i.tex", z, l)

   set key width -3.75

   set title sprintf("Total and Partial ICS Curves [$z_{\\rm{proj}} = %i$]", z)
   set xlabel "Energy [eV]"
   set ylabel "ICS [$\\rm{a}_{0}^{2}$]"

   set xrange [0:50]
   set yrange [0.00001:*]
   set logscale y
   set format y "$10^{%L}$"

   plot \
     ics(z,l) using 1:2 \
       title sprintf("$\\scriptscriptstyle\\sum_{\\ell}$", j-3) \
       with lines lc "black", \
     for [j=3:3+l] \
       ics(z,l) using 1:j \
       title sprintf("$\\scriptscriptstyle\\ell = %i$", j-3) \
       with lines palette frac ((j-3)/(1.0*(l+1)))

   set output
}


# figure: dcs for variety of energies
do for [i=1:n_z] {
   z = z(i)
   l = 5

   set key outside top right width -7

   set output sprintf("figure_dcs_z-%+i_l-%i.tex", z, l)

   set title "DCS Curves"
   set xlabel "$\\theta$ [${\\phantom{.}}^{\\circ}$]"
   set ylabel "DCS [$\\rm{a}_{0}^{2}.\\rm{sr}^{-1}$]"

   set xrange [0:180]
   set yrange [*:*]
   set logscale y
   set format y "%g"

   plot \
     for [j=1:n_en] \
       dcs(z,l, en(j)) using 1:2 \
       title sprintf("$\\scriptscriptstyle %.1f\\phantom{.}\\rm{eV}$", en(j)) \
       with lines palette frac ((j - 1.0)/(n_en - 1.0))

   set output
}
