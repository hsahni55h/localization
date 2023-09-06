/* This robot is located inside of a 1D hallway which has three doors. 
The robot doesn't know where it is located in this hallway, but it has sensors onboard that can tell it, with some amount of precision, whether it is standing in front of a door, or in front of a wall. 
The robot also has the ability to move around - with some precision provided by its odometry data. 
Neither the sensors nor the movement is perfectly accurate, but the robot aims to locate itself in this hallway.

The mobile robot is now moving in the 1D hallway and collecting odometry and perception data. 
With the odometry data, the robot is keeping track of its current position. Whereas, with the perception data, the robot is identifying the presence of doors.

In this quiz, we are aiming to calculate the state of the robot, given its measurements. This is known by the belief: P(Xt|Z)!

Given:

P(POS): The probability of the robot being at the actual position
P(DOOR|POS): The probability of the robot seeing the door given that it’s in the actual position
P(DOOR|¬POS): The probability of the robot seeing the door given that it’s not in the actual position
Compute:

P(POS|DOOR): The belief or the probability of the robot being at the actual position given that it’s seeing the door.*/


#include <iostream>
using namespace std;

int main() {
	
	//Given P(POS), P(DOOR|POS) and P(DOOR|¬POS)
	double a = 0.0002 ; //P(POS) = 0.002
	double b = 0.6    ; //P(DOOR|POS) = 0.6
	double c = 0.05   ; //P(DOOR|¬POS) = 0.05
	
	//TODO: Compute P(¬POS) and P(POS|DOOR)
	double d =   1-a;                //P(¬POS)
	double e = b*a/((a*b) + (d*c));               //P(POS|DOOR)
	
	//Print Result
	cout << "P(POS|DOOR)= " <<    e    << endl;
	
	return 0;
}