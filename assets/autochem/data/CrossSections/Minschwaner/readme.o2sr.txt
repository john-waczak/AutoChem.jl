
In this directory there are three files, each 
contains the coefficients appropriate for the 3 temperatures ranges
as outlined in the 1992 JGR paper.  
"fitcoef_cold" for 130 to 190 K, "fitcoef_mid" for 190 to 280 K, 
and "fitcoef_hot" for 280 to 500 K.  The columns labelled "a0,a1,a2" 
contain the polynomial coefficients such that the cross section = 
1.e-20*(a0*d**2 + a1*d + a2) in units cm**-2, where d is the transformed 
temperature variable defined in our eq. (6):   d = [(T-100)/10]**2.
The last 2 columns list the MAXIMUM percentage error in the cross
section at that wavenumber and the temperature at which this
maximum error occurs.

Please note: The underlying Herzberg continuum must be added
for atmospheric applications.  Also, the poly fits go haywire
if used outside the intended temperature range.

Good luck, and let me know if you're successful in getting the
files.  If you have any questions or problems, feel free to
call (505) 835-5226, or e-mail krm@kestrel.nmt.edu.  

Ken Minschwaner
Department of Physics
New Mexico Institute of Mining and Technology
Socorro, NM
87801
