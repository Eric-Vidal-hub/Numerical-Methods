function barrierPotential = BarrierPotential(positions, barrierStart, barrierEnd, energy)
    % This function computes the barrier potential for a given set of positions.
    % It assigns a constant energy value to positions within the barrier and zero elsewhere.
    %
    % Inputs:
    % positions: A vector of positions
    % barrierStart: The start position of the barrier
    % barrierEnd: The end position of the barrier
    % energy: The energy value for positions within the barrier
    %
    % Output:
    % barrierPotential: A vector of the same size as positions, containing the computed barrier potential for each position

    % Initialize the barrier potential array with zeros
    barrierPotential = zeros(size(positions));
    
    % Assign energy to elements of barrierPotential where corresponding position is within the barrier
    barrierPotential(barrierStart <= positions & positions <= barrierEnd) = energy;
end