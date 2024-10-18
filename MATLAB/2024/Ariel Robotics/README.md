# Aerial Robotics

## Overview
My experience with Aerial Robotics (sadly, the folder was misspelled upon creation) constituted a comprehensive introduction to robotic aircraft. The course examined rotorcraft dynamics modeling, feedback control, sensing, state estimation, path planning, machine vision, and decision-making under uncertainty. Furthermore, we design an automation protocol, written in C++, that commands a quadcopter competing in a game.

## Tournament
The tournament took the form of a capture-the-flag-like game called “pop the balloon,” in a speed trial or a race. Worked in teams of 2 to 3 to develop an optimal control framework.

## Projects
To work towards our final competition, each project built on complexity and computational demands.
In the earlier labs, we individually developed a simulator for a UAV, including simulating the UAV dynamics and the sensor measurement. We then processed simulated measurements to estimate the UAV attitude and position. Exploiting the position and attitude estimates, each closed the attitude and position control loops within their own simulator. At the end, we implemented some of our code in real time within the feedback loop of a live quad-rotor UAV. In later projects, each member also developed mapping and path planning methods.
