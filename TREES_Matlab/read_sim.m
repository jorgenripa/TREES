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
    load([simpath resultsName '.mat'],'sim');
    return
end
if isempty(fsim) % there is no file
    sim = NaN;
    return
end

fid = fopen([simpath resultsName '.sim'],'r');

file_version = fread(fid,1,'int32');
if ~ismember(file_version,[1 2 3 4])
    fclose(fid);
    error(['Wrong file version: ' num2str(file_version)])
end

if file_version < 3
    sim.seed = fread(fid,1,'uint32=>uint32');
else
    sim.seed = fread(fid,1,'uint64=>uint64');
end

% Read parameters

% Simulation:
if file_version < 3
    verbose = readPar(fid,'verbose');
end
sim.t_max = readPar(fid,'t_max');
sim.sample_interval = readPar(fid,'sample_interval');
if file_version >= 3
    sim.microsamples_option = readPar(fid,'microsamples');
end
if file_version >= 2
    sim.checkpoint_interval = readPar(fid,'checkpoint_interval');
    [parname,parval] = readParPair(fid);
    if strcmp(parname,'keep_old_checkpoints')
        sim.keep_old_checkpoints = parval;
        seed2 = readPar(fid,'seed');
    elseif strcmp(parname,'seed')
        if isnan(str2double(parval))
            seed2 = parval;
        else
            seed2 = str2double(parval);
        end
    else
        error(['Expected : keep_old_checkpoints, found : ' parname])
    end
else
    seed2 = readPar(fid,'seed');
end
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
if strcmpi(sim.Genetics.model,'omnigenic')
    sim.Genetics.loci = readPar(fid,'loci');
end

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
    sim.Fitness{end+1} = readFitness(modName,fid,file_version);
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
if file_version >= 3
    sim.microsamples = read_microsamples(fid, sim, file_version );
end
sim.sample_count = fread(fid,1,'int32');
sample0 = readSample(fid,sim,file_version);
sim.samples(1) = sample0;
sim.samples(sim.sample_count) = sample0; % preallocation
for si=2:sim.sample_count
    sim.samples(si) = readSample(fid,sim,file_version);
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
        next_id = fread(fid,1,'ubit64'); % This is not used
        sim.gene_lists{li} = list;
    end
end
% Don't read checkpoints
fclose(fid);
% Make simple stats:
sim.stats = calc_stats(sim);
save([simpath resultsName '.mat'], 'sim','-v7.3')


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
    case 'normal_deviate'
        transform.SD = readPar(fid,'sd');
    case 'range'
        transform.min = readPar(fid,'min');
        transform.max = readPar(fid,'max');
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
        if lower(space.dispersal_type(1)) == 'd'
            space.dispersal_distance = readPar(fid,'dispersal_distance');
        end
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

function fitness = readFitness(module,fid, file_version)
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
        if file_version >= 2
            fitness.optimal_value = readPar(fid,'optimal_value');
        end
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
    case 'catastrophes'
        fitness.P_catastrophe = readPar(fid,'p_catastrophe');
        fitness.P_survive = readPar(fid,'p_survive');
    otherwise
        error('Unknown Fitness module!')
end

function ms = read_microsamples(fid, sim,file_version)
if file_version==3
    microsamples_option = fread(fid,1,'char');
else
    microsamples_option = [];
end
microsamples_count = fread(fid,1,'int32');
if microsamples_count == 0
    ms = [];
else
    tot_dim = 0;
    for ti=1:length(sim.Traits)
        if sim.Traits(ti).loci_per_dim >0
            tot_dim = tot_dim + sim.Traits(ti).dims;
        end
    end
    if ~strcmpi(sim.Space.model,'none')
        tot_dim = tot_dim + sim.Space.dimensions;
    end
    % pre-allocate matrices:
    T = zeros(1, microsamples_count,'int64');
    CPU_TIME = zeros(1, microsamples_count,'single');
    N = zeros(1, microsamples_count,'int32');
    MEANS = zeros(tot_dim, microsamples_count,'single');
    VARS = zeros(tot_dim, microsamples_count,'single');
    COVARS = zeros(tot_dim*(tot_dim+1)/2, microsamples_count,'single');
    
    for msi=1:microsamples_count
        if file_version > 3
            microsamples_option = fread(fid,1,'char=>char');
        end
        T(msi) = fread(fid,1,'int64');
        CPU_TIME(msi) = fread(fid,1,'double');
        N(msi) = fread(fid,1,'int32');
        switch microsamples_option
            case 'm'
                if file_version>3
                    msize = fread(fid,1,'int32');
                end
                MEANS(:,msi) = fread(fid,tot_dim,'float32=>float32');
                if file_version>3
                    vsize = fread(fid,1,'int32'); % This should be zero
                end
            case 'v'
                if file_version>3
                    msize = fread(fid,1,'int32');
                end
                MEANS(:,msi) = fread(fid,tot_dim,'float32=>float32');
                if file_version>3
                    vsize = fread(fid,1,'int32');
                end
                VARS(:,msi) = fread(fid,tot_dim,'float32=>float32');
            case 'c'
                if file_version>3
                    msize = fread(fid,1,'int32');
                end
                MEANS(:,msi) = fread(fid,tot_dim,'float32=>float32');
                if file_version>3
                    vsize = fread(fid,1,'int32');
                end
                COVARS(:,msi) = fread(fid,tot_dim*(tot_dim+1)/2,'float32=>float32');
        end
    end
    ms.T = T;
    ms.N = N;
    ms.CPU_TIME = CPU_TIME;    
    ms.MEANS = MEANS;
    if any(VARS(:)>0)
        ms.VARS = VARS;
    end
    if any(COVARS(:)~=0)
        ms.COVARS = COVARS;
    end
end

function sample = readSample(fid,sim, file_version)
sample.gen = fread(fid,1,'int64');
if file_version >= 3
    sample.cputime = fread(fid,1,'float32');
end
sample.size = fread(fid,1,'int32'); % file format is little endian.
if file_version > 3
    gene_sampling = fread(fid,1,'int8'); %#ok<NASGU> % extra stored bool
end
if sim.gene_sampling
    if file_version > 3
        gene_type = fread(fid,1,'int8'); %#ok<NASGU> % extra stored type char
    end
    if file_version > 3 
        loci = fread(fid,1,'int32');
        pop_size = fread(fid,1,'int32');
        gene_tracking = fread(fid,1,'int8');
        if sim.gene_tracking
            m_n = fread(fid,2,'int32');
            G1id = fread(fid,sim.loci*sample.size,'ubit64=>ubit64');
            m_n = fread(fid,2,'int32');
            G2id = fread(fid,sim.loci*sample.size,'ubit64=>ubit64');
        end
    end
    switch lower(sim.Genetics.model)
        case 'diallelic'
            if file_version > 3
                G1 = read_bit_vector(fid);
                G2 = read_bit_vector(fid);
            else
                G1 = fread(fid,sim.loci*sample.size,'int8=>int8');
                G2 = fread(fid,sim.loci*sample.size,'int8=>int8');
            end
        case {'continuous_alleles', 'omnigenic'}
            if file_version > 3, m_n = fread(fid,2,'int32'); end
            G1 = fread(fid,sim.loci*sample.size,'float32=>float32');
            if file_version > 3, m_n = fread(fid,2,'int32'); end
            G2 = fread(fid,sim.loci*sample.size,'float32=>float32');
    end
    sample.G1 = reshape(G1,sim.loci,sample.size);
    sample.G2 = reshape(G2,sim.loci,sample.size);
end
if file_version>=2
    if sim.gene_tracking
        if file_version==3
            G1id = fread(fid,sim.loci*sample.size,'ubit64=>ubit64');
            G2id = fread(fid,sim.loci*sample.size,'ubit64=>ubit64');
        end
        sample.G1id = reshape(G1id,sim.loci,sample.size);
        sample.G2id = reshape(G2id,sim.loci,sample.size);
    end
end

if file_version>3
    ntraits = fread(fid,1,'int32');
end
for tri = 1:length(sim.Traits)
    tr = sim.Traits(tri);
    if tr.loci_per_dim > 0
        if file_version>3, m_n = fread(fid,2,'int32'); end
        sample.(tr.name) = reshape( fread(fid,sample.size*tr.dims,'float32'), tr.dims, sample.size);
    end
end

% Space:
if file_version>3
    space_type = fread(fid,1,'int8');
    size_dims = fread(fid,2,'int32');
    tot_length = fread(fid,1,'int32'); % should be prod(size_dims)
end
switch lower(sim.Space.model)
    case 'none'
    case 'discrete' 
        L = sim.Space.size;
        D = sim.Space.dimensions;
        if file_version>3 % linear_patches stored
            linear_patches = fread(fid,sample.size,'int32');
            sample.pos = zeros(D, sample.size);
            for d=D:-1:1
                sample.pos(d,:) = mod(linear_patches, L);
                linear_patches = floor(linear_patches/L);
            end
        else
            sample.pos = reshape(fread(fid,sample.size*D,'int32'),D,sample.size);
        end
    case 'continuous'
        sample.pos = reshape(fread(fid,sample.size*sim.Space.dimensions,'float32'),sim.Space.dimensions,sample.size);
    otherwise
        error('Unsupported space model')
end
if file_version==1
    if sim.gene_tracking
        G1id = fread(fid,sim.loci*sample.size,'ubit64');
        sample.G1id = reshape(G1id,sim.loci,sample.size);
        G2id = fread(fid,sim.loci*sample.size,'ubit64');
        sample.G2id = reshape(G2id,sim.loci,sample.size);
    end
end

function st = calc_stats(sim)
ndims = 0;
names = {};
for tr=sim.Traits
    if tr.loci_per_dim>0
        ndims = ndims + tr.dims;
        if tr.dims==1
            names{end+1} = tr.name;
        else
            for d=1:tr.dims
                names{end+1} = [tr.name '_' num2str(d)];
            end
        end
    end
end
if ~strcmp(sim.Space.model,'none')
    ndims = ndims+sim.Space.dimensions;
    names{end+1} = 'pos';
end

MM = zeros(ndims,sim.sample_count);
CC = zeros(ndims*(ndims+1)/2, sim.sample_count);
for si=1:sim.sample_count
    sa = sim.samples(si);
    X = [];
    for tr=sim.Traits
        if tr.loci_per_dim>0
            X = [X sa.(tr.name)'];
        end
    end
    if ~strcmp(sim.Space.model,'none')
        X = [X sa.pos'];
    end
    MM(:,si) = mean(X)';
    C = cov(X,1); % complete population sample
    CC(:,si) = C(tril(ones(size(C),'logical')));
end
st.names = names;
st.means = MM;
st.covars = CC;

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

function v = read_bit_vector(fid)
bit_count = fread(fid,1,'int32');
chunk_size = fread(fid,1,'int32'); % should be 8 (size in bytes)
chunk_bits = chunk_size*8; % should be 64
chunk_count = ceil(bit_count/chunk_bits);
chunks = fread(fid,chunk_count,'ubit64=>ubit64');
if exist('bits2bytes.mexmaci64','file')
    v = bits2bytes(chunks,bit_count);
else
    v = zeros(1,bit_count,'int8');
    pos = 1;
    for ch_i = 1:floor(bit_count/chunk_bits)
        chunk = chunks(ch_i);
        v(pos:(pos+63)) = bitget(chunk,1:64);
        pos = pos+64;
        % This is slower:
        %     mask = uint64(1);
        %     for bit=1:chunk_bits
        %         v(pos) = bitand(mask,chunk)>0;
        %         pos = pos+1;
        %         mask = mask*uint64(2);
        %     end
    end
    bits_left = rem(bit_count,chunk_bits);
    if bits_left>0
        chunk = chunks(end);
        v(pos:(pos+bits_left-1)) = bitget(chunk,1:bits_left);
    end
end
v = v*2 - 1; % convert (0/1) to (-1/+1)
