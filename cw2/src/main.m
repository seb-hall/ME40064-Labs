%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ME40064 Coursework 2
%
% File         :  main.m
% Author       :  11973
% Created      :  2025-11-24 (YYYY-MM-DD)
% License      :  MIT
% Description  :  Main function for solving transient diffusion equation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main()
    fprintf("ME40064 Coursework 2 Starting...\n");
    
    % Part 1: Software Verification
    Coursework.Part1Plots();
    Coursework.Part1Convergence();

    % Part 2: Software Features
    Coursework.Part2TimeIntegrationComparison();
    Coursework.Part2GaussianQuadrature();

    % Part 3: Modelling & Simulation Results
    Coursework.Part3InitialResults();
    Coursework.Part3MinimumEffectiveDose();
    Coursework.Part3DoseSensitivityAnalysis();

    fprintf("...ME40064 Coursework 2 Complete\n");
end
