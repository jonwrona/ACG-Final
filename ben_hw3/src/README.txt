HOMEWORK 3: RAY TRACING, RADIOSITY, & PHOTON MAPPING

NAME:  < Benjamin West >


TOTAL TIME SPENT:  < Umm no clue probably between 10 and low 20s... >
Please estimate the number of hours you spent on this assignment.


COLLABORATORS: 
You must do this assignment on your own, as described in the 
Academic Integrity Policy.  If you did discuss the problem or errors
messages, etc. with anyone, please list their names here.



RAYTRACING:
< insert comments on the implementation, known bugs, extra credit >
I've got raytracing fully functional, and for soft shadows I implemented both random selection of points and always choosing the 4 same points. I personally prefer the second. I made it very easy to switch between the two implementations with comments on raytracer.cpp line 125



RADIOSITY:
< insert comments on the implementation, known bugs, extra credit >
I used the Southwell relaxation solver. I also added code to distribute the undistributed light ambiently as it was solved. This slightly messed up my undistributed light visualization (which makes sense...). Anyway if you want to comment it in you can on line
254 of radiosity.cpp



PHOTON MAPPING:
< insert comments on the implementation, known bugs, extra credit >
As far as I know this is fully functional with no known bugs. Takes a while to render, but that is to be expected



OTHER NEW FEATURES OR EXTENSIONS FOR EXTRA CREDIT:
Include instructions for use and test cases and sample output as appropriate.
