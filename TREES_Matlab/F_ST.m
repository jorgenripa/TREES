function [F_ST_TOT, F_ST_i] = F_ST(sample, group1, group2, loci)

loci = loci(:)'; % make sure it's a row vector

r = 2; % number of groups

if ~isfield(sample,'G1')
    error('Gene sampling has to be turned on to calculate F_ST')
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
n_bar = (n1+n2)/2;
n_c = (N - (n1^2+n2^2)/N);

L = length(loci);
F_ST_i = zeros(L,1);
% Following Weir & Cockerham (1984):
a_tot = 0;
b_tot = 0;
c_tot = 0;

for li = 1:length(loci)
    locus = loci(li);
    % group genes:
    G1 = [sample.G1(locus,group1); sample.G2(group1)];
    G2 = [sample.G1(locus,group2); sample.G2(group2)];
    a_locus = 0;
    b_locus = 0;
    c_locus = 0;
    
    % Find all alleles:
    all_alleles = unique([G1(:);G2(:)])';
    for allele = all_alleles
        p1 = sum(G1(:)==allele)/2/n1;
        p2 = sum(G2(:)==allele)/2/n2;
        p_bar = (n1*p1+n2*p2)/N;
        s2 = (n1*(p1-p_bar)^2 + n2*(p2-p_bar)^2)/(r-1)/n_bar; 
        % homozygotes:
        homo1 = sum(G1(1,:)==allele & G1(2,:)==allele);
        % proportion heterozygotes:
        hetero1 = (p1*2*n1 - 2*homo1) / n1;
        % homozygotes:
        homo2 = sum(G2(1,:)==allele & G2(2,:)==allele);
        % proportion heterozygotes:
        hetero2 = (p2*2*n2 - 2*homo2) / n2;
        h_bar = (n1*hetero1 + n2*hetero2)/N; % average heterozygosity
        a = n_bar/n_c*( s2 - 1/(n_bar-1)*( p_bar*(1-p_bar) - (r-1)/r*s2 - h_bar/4));
        a_locus = a_locus + a;
        b = n_bar/(n_bar-1)*( p_bar*(1-p_bar) - (r-1)/r*s2 - (2*n_bar-1)/4/n_bar * h_bar );
        b_locus = b_locus + b;
        c = h_bar/2;
        c_locus = c_locus + c;
    end
    F_ST_i(li) = a_locus/(a_locus + b_locus + c_locus);
    a_tot = a_tot + a_locus;
    b_tot = b_tot + b_locus;
    c_tot = c_tot + c_locus;
end
F_ST_TOT = a_tot/(a_tot + b_tot + c_tot);
