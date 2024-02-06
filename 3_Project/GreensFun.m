function [greenFunction]= GreensFun(energyStep,position1,position2,energy,damping,energyScale)
    %% GREEN'S FUNCTION
    % Function to calculate the Green's function for given parameters and positions.
    %
    % Syntax:  greenFunction = Greenf(energyStep, position1, position2, energy, damping, energyScale)
    %
    % Inputs:
    %    energyStep - A parameter related to the energy level
    %    position1, position2 - The positions for which to calculate the Green's function
    %    energy - The energy level
    %    damping - A damping factor
    %    energyScale - A scaling factor for the energy
    %
    % Outputs:
    %    greenFunction - Calculated Green's function

    % Calculate the wave number outside the potential
    waveNumberOutside = sqrt((energy + 1i * damping) / energyScale);
    
    % Calculate the wave number inside the potential
    waveNumberInside = sqrt((energy + 1i * damping - energyStep) / energyScale);
    
    % Calculate the common term used in the Green's function
    commonTerm = 1 / (2i * waveNumberInside);

    % If both positions are greater than or equal to 0, calculate the Green's function using equation [A.76]
    if(position1 >= 0 && position2 >= 0)
        greenFunction = commonTerm * (exp(1i * waveNumberInside * abs(position1 - position2)) +  exp(1i * waveNumberInside * (position1 + position2)) * ((waveNumberInside - waveNumberOutside) / (waveNumberInside + waveNumberOutside)));
    end

    % If position1 and position2 are on different sides of the potential, calculate the Green's function using equation [A.77]
    if (position1 < 0 && position2 >= 0) || (position2 < 0 && position1 >= 0)
        greenFunction = exp(-1i * waveNumberOutside * min(position1, position2) + 1i * waveNumberInside * max(position1, position2)) / (1i * (waveNumberOutside + waveNumberInside));
    end
end