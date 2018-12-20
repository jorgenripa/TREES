function sim = read_sim(parfilepath,index)

pathpos = strfind(parfilepath,'/');
if ~isempty(pathpos)
    pathpos = pathpos(end);
else
    pathpos = 0;
end
simpath = parfilepath(1:pathpos);
sim.name = parfilepath(pathpos+1:end-4); % remove extension

resultsName = [sim.name '_results_' num2str(index)];% '.sim'];
fmat = dir([simpath resultsName '.mat']);
fsim = dir([simpath resultsName '.sim']);

if ~isempty(fmat) && (isempty(fsim) || fsim.datenum<fmat.datenum)
    S = load([simpath resultsName '.mat']);
    sim = S.sim;
    return
end
if isempty(fsim) % there is no file
    sim = NaN;
    return
end

fid = fopen([simpath resultsName '.sim'],'r');

fileVersion = fread(fid,1,'int32');
sim.seed = fread(fid,1,'uint32');
if fileVersion~=1
    fclose(fid);
    if fileVersion == 6
        sim = readSim7(parfilepath,index);
        return
    else
        disp(fileVersion)
        error('Wrong file version:')
    end
end

% Read parameters

% Simulation:
verbose = readPar(fid,'verbose');
sim.t_max = readPar(fid,'t_max');
sim.sample_interval = readPar(fid,'sample_interval');
seed2 = readPar(fid,'seed');
%if seed2 ~= sim.seed
%    warning('Seed error')
%end
sim.gene_tracking = readPar(fid,'gene_tracking')=='Y';
sim.gene_sampling = readPar(fid,'gene_sampling')=='Y';

%Population
sim.F = readPar(fid,'f');
sim.n_0 = readPar(fid,'n_0');

%Genetics
sim.Genetics.model = readPar(fid,'genetics');
sim.Genetics.P_mutation = readPar(fid,'p_mutation');

%Traits:
sim.Traits = struct('name',{},'dims',{},'loci_per_dim',{},'initial_value',{},'transforms',{});
[s1,s2] = readParPair(fid);
while strncmp(s1,'trait',5)
    trait = [];
    trait.name = s2;
    trait.dims = readPar(fid,'dimensions');
    if strcmp(s1,'trait')
        trait.loci_per_dim = readPar(fid,'loci_per_dim');
    else
        trait.loci_per_dim = 0; % it's a constant
    end
    trait.initial_value = readPar(fid,'initial_value');
    trait.transforms = {};
    [s1,s2] = readParPair(fid);
    while strcmp(s1,'transform')
        trait.transforms{end+1} = readTransform(s2,fid);
        [s1,s2] = readParPair(fid);
    end
    sim.Traits(end+1) = trait;
end

modType = s1;
modName = s2;
% read Space module if found:
if strcmp(modType,'space')
    sim.Space = readSpace(modName,fid);
    [modType,modName] = readParPair(fid);
else
    sim.Space.model = 'none';
end
%read fitness modules:
sim.Fitness = {};
while strcmp(modType,'fitness')
    sim.Fitness{end+1} = readFitness(modName,fid);
    [modType,modName] = readParPair(fid);
end

if strcmp(modType,'mating_pool')
    sim.Mating.pool = modName;
    if lower(modName(1))=='l' % local mating
        sim.Mating.s_space = readPar(fid,'s_space');
    end
else
    error(['Mating_pool :' modType ' : ' modName])
end
sim.Mating.trials = readPar(fid,'mating_trials');
sim.Mating.preferences = {};
[modType,modName] = readParPair(fid);
while strcmp(modType,'mating_preference')
    sim.Mating.preferences{end+1} = readPreference(modName,fid);
    [modType,modName] = readParPair(fid);
end

%%% end parameter file

sim.loci = fread(fid,1,'int32');
sim.sample_count = fread(fid,1,'int32');
for si=1:sim.sample_count
    sim.samples(si) = readSample(fid,sim);
end
if sim.gene_tracking
    sim.gene_lists = {};
    for li=1:sim.loci
        listsize = fread(fid,1,'int32');
        list = struct('id',cell(1,listsize),'parent',0,'birth',0,'death',0,'effect',0,'children',0,'child_list',0);
        for ai=1:listsize
            a = [];
            a.id = fread(fid,1,'ubit64');
            a.parent = fread(fid,1,'ubit64');
            a.birth = fread(fid,1,'int64');
            a.death = fread(fid,1,'int64');
            a.effect = fread(fid,1,'float32');
            a.children = fread(fid,1,'int32');
            a.child_list = fread(fid,a.children,'ubit64');
            list(ai) = a;
        end
        sim.gene_lists{li} = list;
    end
end

fclose(fid);
save([simpath resultsName '.mat'], 'sim')


function transform = readTransform(name,fid)
transform.name = name;
switch lower(name)
    case 'linear'
        transform.offset = readPar(fid,'offset');
        transform.scale = readPar(fid,'scale');
    case 'abs'
        % no parameters
    case 'logistic'
        transform.min = readPar(fid,'min');
        transform.max = readPar(fid,'max');
    case 'normal_deviation'
        transform.SD = readPar(fid,'sd');
    otherwise
        error(['Unknown transform : ' name])
end

function space = readSpace(module,fid)
space.model = module;
switch lower(module)
    case 'none'
        % No parameters
    case 'discrete'
        space.size = readPar(fid,'size');
        space.dimensions = readPar(fid,'dimensions');
        space.P_disperse = readPar(fid,'p_disperse');
        space.dispersal_type = readPar(fid,'dispersal_type');
        space.boundary = readPar(fid,'boundary');
        space.initial_position = readPar(fid,'initial_position');
    case 'continuous'
        space.size = readPar(fid,'size');
        space.dimensions = readPar(fid,'dimensions');
        space.P_disperse = readPar(fid,'p_disperse');
        space.dispersal_distance = readPar(fid,'dispersal_distance');
        space.boundary = readPar(fid,'boundary');
        space.initial_position = readPar(fid,'initial_position');
    otherwise
        error(['Unknown Space model : ' module])
end

function mating = readPreference(module,fid)
mating.name = module;
switch lower(module)
    case 'target_selection'
        mating.display = readPar(fid,'display');
        mating.preference = readPar(fid,'preference');
        mating.strength = readPar(fid,'strength');
        mating.disassortative_limit = readPar(fid,'disassortative_limit');
    otherwise
        error('Unknown Preference model!')
end

function fitness = readFitness(module,fid)
fitness.name = module;
switch lower(module)
    case 'resource_landscape'
        fitness.trait = readPar(fid,'trait');
        fitness.r = readPar(fid,'r');
        fitness.K_0 = readPar(fid,'k_0');
        fitness.s_K = readPar(fid,'s_k');
        fitness.s_a = readPar(fid,'s_a');
        fitness.s_space = readPar(fid,'s_space');
        fitness.k_space = readPar(fid,'k_space');
    case 'stabilizing_selection'
        fitness.trait = readPar(fid,'trait');
        fitness.cost_coefficient = readPar(fid,'cost_coefficient');
        fitness.cost_exponent = readPar(fid,'cost_exponent');
    case {'discrete_resources'}
        fitness.trait = readPar(fid,'trait');
        fitness.n_R = readPar(fid,'n_r');
        fitness.K = readPar(fid,'k');
        fitness.a_0 = readPar(fid,'a_0');
        fitness.s_a = readPar(fid,'s_a');
        fitness.c_min = readPar(fid,'c_min');
    case 'spatial_gradient'
        fitness.trait = readPar(fid,'trait');
        fitness.k_space = readPar(fid,'k_space');
        fitness.s_selection = readPar(fid,'s_selection');
    case 'density_dependence'
        fitness.r = readPar(fid,'r');
        fitness.K = readPar(fid,'k');
        fitness.s_space = readPar(fid,'s_space');
    case 'catastrophe'
        fitness.P_catastrophe = readPar(fid,'p_catastrophe');
        fitness.P_survive = readPar(fid,'p_survive');
    otherwise
        error('Unknown Fitness module!')
end


function sample = readSample(fid,sim)
sample.gen = fread(fid,1,'int64');
sample.sample_size = fread(fid,1,'int32'); % file format is little endian.

if sim.gene_sampling
    switch lower(sim.Genetics.model)
        case 'diallelic'
            G1 = fread(fid,sim.loci*sample.sample_size,'int32');
            G2 = fread(fid,sim.loci*sample.sample_size,'int32');
        case 'continuous_alleles'
            G1 = fread(fid,sim.loci*sample.sample_size,'float32');
            G2 = fread(fid,sim.loci*sample.sample_size,'float32');
    end
    sample.G1 = reshape(G1,sim.loci,sample.sample_size);
    sample.G2 = reshape(G2,sim.loci,sample.sample_size);
end

for tri = 1:length(sim.Traits)
    tr = sim.Traits(tri);
    if tr.loci_per_dim > 0
        sample.(tr.name) = reshape( fread(fid,sample.sample_size*tr.dims,'float32'), tr.dims, sample.sample_size);
    end
end

switch lower(sim.Space.model)
    case 'none'
    case 'discrete'
        sample.pos = reshape(fread(fid,sample.sample_size*sim.Space.dimensions,'int32')',sim.Space.dimensions,sample.sample_size);
    case 'continuous'
        sample.pos = reshape(fread(fid,sample.sample_size*sim.Space.dimensions,'float32')',sim.Space.dimensions,sample.sample_size);
    otherwise
        error('Unsupported space model')
end
if sim.gene_tracking
    G1id = fread(fid,sim.loci*sample.sample_size,'ubit64');
    sample.G1id = reshape(G1id,sim.loci,sample.sample_size);
    G2id = fread(fid,sim.loci*sample.sample_size,'ubit64');
    sample.G2id = reshape(G2id,sim.loci,sample.sample_size);
end

function val = readPar(fid, pname)
sline = readParLine(fid);
if ~isnan(sline)
    if strncmp(lower(sline),pname,length(pname))
        % find colon:
        pos = find(sline==':');
        val = sline(pos+1:end);
        if ~isnan(str2double(val))
            val = str2double(val);
        end
    else
        error(['Looking for ' pname ', found ' sline])
    end
else
    error(['Looking for ' pname ', found EOF'])
end

function [name,val] = readParPair(fid)
sline = readParLine(fid);
if isnan(sline(1)) % EOF?
    name = '';
    val = NaN;
else
    % find colon:
    pos = find(sline==':');
    name = lower(sline(1:pos-1));
    val = sline(pos+1:end);
end

function s = readParLine(fid)
s = readLine(fid);
while isempty(s)
    s = readLine(fid);
end


function s = readLine(fid)
s1 = fread(fid,1);
if s1==0 % EOF?
    s = NaN;
else
    s = '';
    while ~ismember(s1,[10,13,0,',','#'])
        s = [s s1];
        s1 = fread(fid,1);
    end
    if s1=='#' % skip rest of line
        dummy = fgetl(fid);
    end
    % strip whitespace:
    s = sscanf(s,'%1s');
    % Strip comment
    cpos = find(s=='#',1,'first');
    if ~isempty(cpos)
        s = s(1:cpos-1);
    end
end

