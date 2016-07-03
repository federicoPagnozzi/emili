#rm ./Tests/STD/*

declare -a improvement=("first" "best")
declare -a neighborhood=("twoExc" "ref")
declare -a perturbation=("ref" "twoExc")
declare -a instance=(
		     "/home/antoniofisk/Desktop/Uni/MasterThesis/Project/build/Instance_V_1.11.xml"
		     "/home/antoniofisk/Desktop/Uni/MasterThesis/Project/build/Instance_V_1.11.xml"
		     "/home/antoniofisk/Desktop/Uni/MasterThesis/Project/build/Instance_V_1.11.xml"
		     "/home/antoniofisk/Desktop/Uni/MasterThesis/Project/build/Instance_V_1.11.xml"
		     "/home/antoniofisk/Desktop/Uni/MasterThesis/Project/build/Instance_V_1.11.xml"
		     "/home/antoniofisk/Desktop/Uni/MasterThesis/Project/build/Instance_V_1.11.xml"
		     "/home/antoniofisk/Desktop/Uni/MasterThesis/Project/build/Instance_V_1.11.xml"
		     "/home/antoniofisk/Desktop/Uni/MasterThesis/Project/build/Instance_V_1.11.xml"
		     "/home/antoniofisk/Desktop/Uni/MasterThesis/Project/build/Instance_V_1.11.xml"
		     "/home/antoniofisk/Desktop/Uni/MasterThesis/Project/build/Instance_V_1.11.xml"

		   
		   )


a=1                  # Integer.

for inst in "${instance[@]}"
do
#  ./EMILI $inst IRP IRP ils first random locmin twoExc 1 1 true rndmv twoExc 1 3 1 improve
   ./challengeAL $inst  IRP ils first greedyRandomized 0.1368 1 locmin shiftRem 0.8268 1 0.0733 0.9375 true rndmv ref 0.5727 0.6883 0.6177 1 5 improve -it 1800 rnds $a  
   let "a++"
   
done


#/home/antoniofisk/Desktop/Uni/MasterThesis/Project/build/Instance_V_1.1.xmSl IRP ils first random locmin ref true rndmv twoExc 1 metropolis 3.5 -it 1800



 
