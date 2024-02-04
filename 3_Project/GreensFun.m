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

    waveNumberOutside=sqrt((energy+1i*damping)/energyScale);
    waveNumberInside=sqrt((energy+1i*damping-energyStep)/energyScale);
    commonTerm = 1/(2i*waveNumberInside);

    if(position1>=0 && position2>=0)
        greenFunction=commonTerm * (exp(1i*waveNumberInside*abs(position1-position2)) +  exp(1i*waveNumberInside*(position1+position2))*((waveNumberInside-waveNumberOutside)/(waveNumberInside+waveNumberOutside))); % Green's function [A.76] for both positions greater than 0
    end

    if (position1<0 && position2>=0) || (position2<0 && position1>=0)
        greenFunction=exp(-1i*waveNumberOutside*min(position1,position2)+1i*waveNumberInside*max(position1,position2))/(1i*(waveNumberOutside+waveNumberInside)); % Green's function [A.77] for position1 and position2  on different sides
    end
end