% this function only returns rough values for the dynamic viscosity of air related
% to the temperature given a constant pressure of 1bar (10000Pa) which claims
% that the surrounding fluid (air) pressure of the rocket would not become
% less of its altitude, which it definetly does
% ONLY here due to the absence of dynamic viscosity over the air's pressure
% this rough estimation is used

function ret = dynamicViscosity(temperature)

if(temperature > 25 && temperature < 50){
  ret = 19.67
}
if(temperature > 0 && temperature < 25){
  ret = 18.48
}
if(temperature > -25 && temperature < 0){
  ret = 17.24
}
if(temperature > -50 && temperature < -25){
  ret = 15.96
}
if(temperature > -75 && temperature < -50){
  ret = 14.62
}
if(temperature > -75 && temperature < -50){
  ret = 13.22
}
if(temperature > -100 && temperature < -75){
  ret = 11.77
}
if(temperature > -150 && temperature < -100){
  ret = 8.65
}
