         UV / Vis spectra of atmospheric constituents
                         Version 1

                           Readme



The following items are changed since printing of the booklet:

-  The library reference (ISBN number) has been changed. 
   The correct ISBN number is on the stickers on the booklet´s
   title page.

-  Due to a transfer error the x (lambda) scale of the COFCl spectra 
   on page 5 is missing. 
   The original scale is: lambda = 200 ... 260 nm


New Filter features:

-  To avoid "Division by zero" errors the input check for 
   Lambda[i] < Lambda[i+1] is added. 
   If the check result is false, the (Lambda[i+1]; sigma[i+1] ) 
   value pair is skipped.

-  For better understanding, the mode indicator has been changed 
   from "p0/1 ok" to "Mode 0/1 ok" and the success indicator 
   from "ok" to "File ok"

-  In any case, the filter requests a folder and a file name for 
   file output by opening a standard FILE SAVE window. 
 