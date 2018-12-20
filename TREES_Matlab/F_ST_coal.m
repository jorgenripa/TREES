function [F_ST_tot, F_ST_i, dist_i] = F_ST_coal(sample, sim, group1, group2, loci)

loci = loci(:)'; % make sure it's a row vector

if ~isfield(sample,'G1id')
    error('Gene tracking has to be turned on to calculate F_ST_coal')
end

if islogical(group1)
    group1 = find(group1);
end
if islogical(group2)
    group2 = find(group2);
end
n1 = length(group1);
n2 = length(group2);
N = n1 + n2;
L = length(loci);
F_ST_i = zeros(L,1);
dist_i = zeros(L,1);
for li=1:length(loci)
    locus = loci(li);
    % ids:
    ids1 = [sample.G1id(locus,group1); sample.G2id(locus,group1) ]; 
    ids2 = [sample.G1id(locus,group2); sample.G2id(locus,group2) ];
    ids_all = [ids1 ids2];
    % All alleles present in the population:
    iid = unique(ids_all(:));

    % Phylogenetic distances:
    D = calcAlleleSeparation(iid,sample,sim,locus);
    na = length(iid); % iid is list of unique alleles
    if na>1
        % Relative frequencies:
        Pa1 = zeros(na,1);
        % separation, population 1:
        for ai=1:na
            Pa1(ai) = sum(ids1(:)==iid(ai))/n1/2;
        end
        % Mean separation:
        dist1 = Pa1'*D*Pa1;
        
        % separation, population 2:
        Pa2 = zeros(na,1);
        for ai=1:na
            Pa2(ai) = sum(ids2(:)==iid(ai))/n2/2;
        end
        % Mean separation:
        dist2 = Pa2'*D*Pa2;
        
        % separation, total population:
        Pa = (n1*Pa1 + n2*Pa2)/N;
        % Mean separation:
        dist_all = Pa'*D*Pa;
        dist_mean = (n1*dist1 + n2*dist2)/N;
        F_ST_i(li) = (dist_all - dist_mean)/dist_all;
        dist_i(li) = Pa1'*D*Pa2;
    end
end
F_ST_tot = mean(F_ST_i);

function D = calcAlleleSeparation(iid,sample,sim,li)
% iid : gene id:s
% D : the corresponding distance matrix

na = length(iid);
if na==1
    D = 0;
    return
end
% Phylogenetic distances:
D = zeros(na);
gl = sim.gene_lists{li};
ids = [gl.id];
for a1i = 1:(na-1)
    [pedigree1, tt1] = findPedigree(iid(a1i),gl);
    for a2i = (a1i+1):na
        id2 = iid(a2i);
        separationTime = sample.gen;
        while ~ismember(id2,pedigree1)
            separationTime = gl(find(ids==id2)).birth;
            id2 = gl(find(ids==id2)).parent;
        end
        i = find(id2==pedigree1);
        if i==1 % allele 2 separated from allel 1
            D(a1i,a2i) = sample.gen - separationTime;
        else % both separated from common ancestor
            D(a1i,a2i) = sample.gen - min(separationTime,tt1(i-1));
        end
        D(a2i,a1i) = D(a1i,a2i);
    end
end

function [pedigree, tt] = findPedigree(id,gl)
ids = [gl.id];
pedigree = id;
i = find(ids==id);
tt = gl(i).birth;
while gl(i).parent>0
    pedigree(end+1) = gl(i).parent;
    i = find(ids==pedigree(end));
    tt = [tt gl(i).birth];
end



