% generate polynomial fits from data in simulation files
% simulation files are all found in data/simulation
% these simulations are performed in Python using FEMwell
function out = polyfitFromFile(simFile, polType)
    load(simFile, 'indexTemps', 'indexLambda', 'indexSimResult');
    [meshTemps, meshLambda] = meshgrid(indexTemps, indexLambda);
    vectorTemps = reshape(meshTemps,[],1);
    vectorLambda = reshape(meshLambda,[],1);
    sampleVector = [vectorTemps, vectorLambda];
    nVector = reshape(indexSimResult.(polType)',[],1); % these transposes are really sneaky!
    out = polyfitn(sampleVector, nVector, 4);
end