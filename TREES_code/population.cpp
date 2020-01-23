//
// TREES -
// A TRait-based Eco-Evolutionary Simulation tool
// Copyright (C) 2017  Jörgen Ripa
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Contact: jorgen.ripa@biol.lu.se
//


#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <functional>

#include "population.h"
#include "mating.h"
#include "genotype_phenotype_map.h"
#include "randgen.h"
#include "fastmath.h"

Population::Population(ParameterFile& pf) {
    // Read general parameters:
    gene_tracking = pf.getBool("gene_tracking");
    gene_sampling = pf.getBool("gene_sampling");

    F = pf.getPositiveInt("f");
    
    n0 = (int)pf.getPositiveInt("n_0");
    n = 0;
    age = 0;
    alive.reserve(10000); // Boolean vector, keeping track of who's alive
    alive.assign(0, true); // empty vector
    fitness.reserve(10000);
    fitness.assign(0,1);
    
    // Read genetics model and parameters:
    readGenetics(pf);
    
    // Read traits:
    readTraits(pf);
    
    // read Space module, if specified:
    if (pf.get_next_name()=="space") {
        addSpace(pf);
    } else {
        space = new Null_space(*this);
    }
    
    // Read fitness modules
    fitnessList.clear();
    while (pf.get_next_name()=="fitness") {
        addFitness(pf);
    }
    
    // Mating
    readMating(pf);
    
    
    if (pf.get_next_name().length()>0) {
        std::cout << "Unexpected module at end : " << pf.get_next_name() << "\n";
        exit(0);
    }
    
    // Set efficiency flags:
    withinPatchMating = false;
    if (space->getDims()==0) {
        withinPatchMating = true;
    }
    if (space->isDiscrete() && theMatingType==Local && mate_s_space==0.0) {
        withinPatchMating = true;
    }
}

void Population::readGenetics( ParameterFile& pf) {
    std::string geneticsModel = pf.getStringLower("genetics");
    if ( geneticsModel == "continuous_alleles") {
        genetics = new ContinuousAlleles(*this, pf);
    } else if ( geneticsModel == "diallelic" ) {
        genetics = new Diallelic(*this, pf);
    } else if ( geneticsModel == "omnigenic" ) {
        genetics = new Omnigenic(*this, pf);
    } else {
        std::cout << "Unknown Genetics model!\n";
        exit(0);
    }
}

void Population::readTraits( ParameterFile& pf) {
    traitList.clear();
    while (pf.get_next_name().substr(0,5)=="trait") {
        if (pf.get_next_name()=="trait") {
            addTrait(pf.get_next_value_string(), pf);
        } else if (pf.get_next_name() =="trait_constant") {
            addTraitConstant(pf.get_next_value_string(), pf);
        } else {
            std::cout << "Unknown Trait type:" << pf.get_next_name() << "\n";
            exit(1);
        }
    }
}

void Population::readMating(ParameterFile& pf) {
    if (pf.get_next_name()=="mating_pool") {
        std::string mod_name = pf.get_next_value_string();
        switch (std::tolower(mod_name[0])) {
            case 's':
                theMatingType = Selfing;
                break;
            case 'g':
                theMatingType = Random;
                break;
            case 'l':
                theMatingType = Local;
                mate_s_space = pf.getDouble("s_space");
                break;
            default:
                std::cout << "Unknown Mating_pool type : " << mod_name << "\n";
                exit(0);
                break;
        }
    }else {
        std::cout << "Expected Mating_pool, found : " << pf.get_next_name() << "\n";
        exit(1);
    }

    mating_trials = pf.getPositiveInt("mating_trials");
    
    // Read preference modules
    matingPreferenceList.clear();
    while (pf.get_next_name()=="mating_preference") {
        addPreference(pf);
    }
}

void Population::allAlive() {
    alive.assign(n, true);
    somebodyDied = false;
}

void Population::kill(int individual) {
    alive.at(individual) = false;
    somebodyDied = true;
}

void Population::addTrait(std::string traitName, ParameterFile& pf){
    traitList.push_back(new Trait(traitName, *this, pf));
}
void Population::addTraitConstant(std::string traitName, ParameterFile& pf){
    traitList.push_back(new TraitConstant(traitName, *this, pf));
}

Trait* Population::findTrait(std::string& tname) {
    for (Trait* t :traitList) {
        if (t->getName() == tname) {
            return t;
        }
    }
    std::cout << "Unknown trait: " << tname << '\n';
    exit(0);
}

void Population::addFitness(ParameterFile &pf) {
    std::string mod_name = pf.get_next_value_string_lower();
    if (mod_name=="stabilizing_selection") {
        fitnessList.push_back(new StabilizingSelection(*this, pf));
    } else if (mod_name=="density_dependence") {
        fitnessList.push_back(new DensityDependence(*this, pf));
    } else if (mod_name == "resource_landscape") {
        fitnessList.push_back(new ResourceLandscape(*this, pf));
    } else if (mod_name=="discrete_resources") {
        fitnessList.push_back(new DiscreteResources(*this, pf));
    } else if (mod_name=="spatial_gradient") {
        fitnessList.push_back(new SpatialGradient(*this, pf));
    } else if (mod_name=="catastrophes") {
        fitnessList.push_back(new Catastrophe(*this, pf));
//    } else if (mod_name=="FinchModel") {
//        fitnessList.push_back(new FinchModel(*this, pf));
    } else {
        std::cout << "Unknown Fitness model : " << mod_name << "\n";
        exit(0);
    }
}

void Population::addPreference(ParameterFile &pf) {
    std::string mod_name = pf.get_next_value_string_lower();
    if ( mod_name == "target_selection") {
        matingPreferenceList.push_back(new Preference(*this, pf));
    } else {
        std::cout << "Unknown Mating_preference model : " << mod_name << "\n";
        exit(1);
    }
}

void Population::addSpace(ParameterFile &pf) {
    std::string mod_name = pf.get_next_value_string_lower();
    if ( mod_name == "none") {
        space = new Null_space(*this);
    } else if ( mod_name == "discrete") {
        space = new Discrete_space_imp(*this, pf);
    } else if ( mod_name == "continuous") {
        space = new Continuous_space(*this, pf);
    } else {
        std::cout << "Unknown Space model!\n";
        exit(0);
    }
}

void Population::makeNextGeneration() {
    reproduce(); // mate, produce offspring, replace parents with children
    space->disperse(); // may include dispersal mortality
    // This is to avoid testing for mortality during fitness calculations:
    if (somebodyDied) {
        compactData();
    }
    survive(); // adult survival.
    ++age;
    // Every now and then prune the Allele tables:
    if (gene_tracking && age % 50 == 0) {
        genetics->pruneGeneLists();
    }
}

void Population::reproduce() {
    int nnew = n*F;
    genetics->prepareNewGeneration(nnew);
    space->prepareNewGeneration(nnew);
    std::vector<int> dadList = findDads();
    nnew = 0; // restart counting
    for (int mom=0; mom<n; ++mom) {
        int dad = dadList[mom];
        if (dad>=0) { // there is an accepted mate?
            for (int child=0; child<F; ++child) {
                // recombine from both parents
                genetics->addChild(mom,dad);
                // inherit spatial position (of mom, usually):
                space->addChild(mom, dad);
                ++nnew;
            }
        }
    }
    // Establish next generation:
    n = nnew;
    genetics->nextGeneration();
    space->nextGeneration();
    generatePhenotypes();
    alive.assign(n, true);
    somebodyDied = false;
    fitness.assign(n,1.0);
}

std::vector<int> Population::findDads() {
    std::vector<int> dads(n,-1);
    if (theMatingType == Selfing) {
        for (int mom=0; mom<n; ++mom) {
            dads[mom] = mom;
        }
    } else if (withinPatchMating) {
        Discrete_space& theSpace = dynamic_cast<Discrete_space&>(getSpace());
        std::vector<int> indList;
        for (int patch=0; patch<theSpace.getPatchCount(); ++patch) {
            if (theSpace.popSizeInPatch(patch)>0) {
                int local_n = theSpace.popSizeInPatch(patch);
                // Generate list of patch inhabitants
                indList.clear();
                for (int i=0; i<size(); ++i) {
                    if (theSpace.getLinearPatch(i) == patch) {
                        indList.push_back(i);
                    }
                }
                // Let each Mom choose from the local list
                for (int momi=0; momi<local_n; ++momi) {
                    int mom = indList[momi];
                    int trialCount = 0;
                    while (trialCount<mating_trials) {
                        // First find a possible candidate
                        int picked = (int)(rand1()*local_n);
                        int candidate = indList[picked];
                        // Next, let mom accept or reject according to her preferences
                        double Pchoose = 1;
                        for (Preference* pp : matingPreferenceList) {
                            Pchoose *= pp->getPartnerWeight(mom, candidate);
                        }
                        if (Pchoose==1 || rand1()<Pchoose) {
                            dads[mom] = candidate;// acceptance!
                            break;
                        }
                        ++trialCount; // rejection, try another mate
                    }
                }
            }
        }
    } else if (space->getDims()==0 || theMatingType == Random) {
        // no space:
        for (int mom=0; mom<n; ++mom) {
            int trialCount = 0;
            while (trialCount<mating_trials) {
                // Pick a random candidate
                int candidate = (int)(rand1()*n);
                // Next, let mom accept or reject according to her preferences
                double Pchoose = 1;
                for (Preference* pp : matingPreferenceList) {
                    Pchoose *= pp->getPartnerWeight(mom, candidate);
                }
                if (Pchoose==1 || rand1()<Pchoose) {
                    dads[mom] = candidate; // acceptance!
                    break; // while loop
                }
                ++trialCount; // rejection, try another mate
            }
        }
    } else {
        // General case
        Fastexp fexp;
        double twoss2 = 2 * mate_s_space * mate_s_space;
        if (/* DISABLES CODE */ (false)) {
            for (int mom=0; mom<n; ++mom) {
                std::vector<double> space_weights;
                space_weights.assign(n, 1);
                // Calculate spatial filter weights:
                for (int male=0; male<n; ++male) {
                    double dist2 = getSpace().getDist2(mom, male);
                    space_weights[male] *= fexp.exp(-dist2/twoss2);
                }
                // calculate cumulative sum:
                std::partial_sum(space_weights.begin(), space_weights.end(), space_weights.begin());
                int trialCount = 0;
                while (trialCount<mating_trials) {
                    // First find a possible candidate
                    // It's a random pick, weighted by spatial filters
                    int candidate = weightedChoiceCumSum(space_weights);
                    // Next, let mom accept or reject according to her preferences
                    double Pchoose = 1;
                    for (Preference* pp : matingPreferenceList) {
                        Pchoose *= pp->getPartnerWeight(mom, candidate);
                    }
                    if (Pchoose==1 || rand1()<Pchoose) {
                        dads[mom] = candidate; // acceptance!
                        break; // while loop
                    }
                    ++trialCount; // rejection, try another mate
                }
            }
        }else { // other solution, considerably faster
            for (int mom=0; mom<n; ++mom) {
                int trialCount = 0;
                while (trialCount<mating_trials) {
                    // First find a possible candidate
                    // It's a random pick, weighted by spatial filters
                    int candidate = -1;
                    while (candidate<0) {
                        int male = int(rand1()*n);
                        double dist2 = getSpace().getDist2(mom, male);
                        if(rand1() < fexp.exp(-dist2/twoss2)) {
                            candidate = male;
                        }
                        
                    }
                    // Next, let mom accept or reject according to her preferences
                    double Pchoose = 1;
                    for (Preference* pp : matingPreferenceList) {
                        Pchoose *= pp->getPartnerWeight(mom, candidate);
                    }
                    if (Pchoose==1 || rand1()<Pchoose) {
                        dads[mom] = candidate; // acceptance!
                        break; // while loop
                    }
                    ++trialCount; // rejection, try another mate
                }
            }
        }
    }
    return dads;
}


void Population::generatePhenotypes() {
    for (Trait* tr : traitList) {
        tr->generatePhenotypes();
    }
    fitness.assign(n, 1);
}

void Population::survive() {
    for (Fitness* sp : fitnessList) {
        sp->aggregateFitness(fitness);
    }
    for (int i=0; i<n; ++i) {
        alive.at(i) = alive.at(i) && rand1()<fitness.at(i)/F;
    }
    compactData();
}

void Population::compactData(){
    // compact all data:
    genetics->compactData(alive);
    for( Trait* tr : traitList) {
        tr->compactData(alive);
    }
    space->compactData(alive);
    // compact fitness:
    int iw = 0;
    for (int ir=0; ir<alive.size(); ++ir) {
        if (alive[ir]) {
            fitness[iw]=fitness[ir];
            ++iw;
        }
    }
    n = iw;
    fitness.resize(n);
    alive.assign(n, true);
    somebodyDied = false;
}

Population::~Population() {
    //std::cout << "killing population\n";
    for( Trait* tr : traitList) {
        delete tr;
    }
    for( Fitness* sp : fitnessList) {
        delete sp;
    }
    delete space;
    for( Preference* pp : matingPreferenceList) {
        delete pp;
    }
    delete genetics;
}

void Population::initialize() {
    n = n0;
    // initialize genomes:
    genetics->initialize();// genome pre-allocation
    
    for( Trait* tr : traitList) {
        tr->initialize();// Normally sets genes to match initial trait values
    }
    if (gene_tracking) {
        // create gene tables and assign initial values:
        genetics->initializeTracking();
    }
    // assign initial position:
    space->initialize();
    alive.assign(n, true);
    generatePhenotypes();
    age = 0;
}

int Population::calc_total_data_dimensions() {
    int totd = 0;
    for (Trait* tr : traitList) {
        if (!tr->is_constant()) {
            totd += tr->get_dims();
        }
    }
    totd += space->getDims();
    return totd;
}

void Population::resumeAtCheckpoint(Population_checkpoint& cp) {
    n = cp.get_pop_size();
    age = cp.get_generation();
    genetics->resumeFromCheckpoint(cp.genetics, cp.geneListsCopy);
    space->resumeFromCheckpoint(cp.space);
    if(cp.GP_map_data.size()>0) {
        int gp_map_index = 0;
        // This only restores parameters (if any).
        // Actual phenotypes are regenerated directly from genotypes (below).
        for (Trait* t : traitList) {
            if (!t->is_constant()) {
                t->get_GP_map()->resume_from_checkpoint_data(cp.get_map_data(gp_map_index));
                gp_map_index++;
            }
        }
    }
    generatePhenotypes();
    alive.assign(n, true);
    somebodyDied = false;
    fitness.assign(n, 1.0);
}

Sample_base::~Sample_base(){
}

/*class Population_sample {
 protected:
 time_type generation;
 double cputime; // in seconds
 int pop_size;
 Genetics_sample* genes_sample;
 std::vector<Trait_sample> traits;
 Space_sample* space_sample;
 public:
 Population_sample(Population& pop, double cputime);
 void write_to_file(oSimfile& osf);
 };
*/

Population_sample::Population_sample(Population& pop, double cputime) :
    Sample_base(pop.getAge(), cputime, pop.size()) {
        
    if (pop.gene_sampling) {
        genes_sample = pop.getGenetics().make_sample();
    } else {
        genes_sample = NULL;
    }
    traits.clear();
    for(Trait* tr : pop.traitList) {
        if (!tr->is_constant()) {
            traits.push_back(Trait_sample(*tr));
        }
    }
    space_sample = pop.getSpace().make_sample();
}

Population_sample::~Population_sample() {
    if (genes_sample) {
        delete genes_sample;
    }
    delete space_sample;
}

void Population_sample::write_to_file(oSimfile &osf) {
    osf.write<time_type>(generation);
    osf.write<float>(cputime);
    osf.write<size_type>(pop_size);

    if(genes_sample) {
        osf.write<char>((char)true); // sizeof(bool) is implementation-specific, but sizeof(char) is always 1
        genes_sample->write_to_file(osf);
    } else {
        osf.write<char>((char)false);
    }

    osf.write<int>((int)traits.size()); // The number of sampled traits (not constants)
    for(Trait_sample& trs : traits) {
        trs.write_to_file(osf);
    }

    space_sample->write_to_file(osf);
}

Population_sample::Population_sample(iSimfile& isf) {
    generation = isf.read<time_type>();
    cputime = isf.read<float>();
    pop_size = isf.read<size_type>();
    bool gene_sampling = (bool)isf.read<char>();
    if (gene_sampling) {
        genes_sample = Genetics_sample::create_from_file(isf);
    } else {
        genes_sample = NULL; // Important!!!
    }
    int ntraits = isf.read<int>();
    traits.clear();
    for (int ti=0; ti<ntraits; ++ti) {
        traits.push_back(Trait_sample(isf));
    }
    space_sample = Space_sample::create_from_file(isf);
}

Population_checkpoint::Population_checkpoint(Population& pop, seed_type cpseed, double cputime) :
Sample_base(pop.getAge(), cputime, pop.size()) {
    seed = cpseed;
    // copy all genomes, and possibly ID:s:
    genetics = pop.getGenetics().make_sample();
    // Store gene lists, if applicable:
    if (pop.isTrackingGenes()) {
        geneListsCopy = pop.getGenetics().getGeneLists(); // deep copy
    } else {
        geneListsCopy.clear();
    }
    GP_map_data.clear();
    for(Trait* tr : pop.traitList) {
        if (!tr->is_constant()) {
            std::vector<double>* data = tr->get_GP_map()->make_checkpoint_data();
            if (data) {
                GP_map_data.push_back(data);
            }
        }
    }
    space = pop.getSpace().make_sample();
}

Population_checkpoint::~Population_checkpoint() {
    delete genetics;
    delete space;
}

void Population_checkpoint::write_to_file(oSimfile &osf) {
    osf.write<seed_type>(seed);
    osf.write<time_type>(generation);
    osf.write<double>(cputime);
    osf.write<int>(pop_size);
    genetics->write_to_file(osf);
    osf.write<size_type>((size_type)geneListsCopy.size());
    //One genelist per locus:
    for (GeneList& list : geneListsCopy) {
        list.writeGenes(osf);
    }
    osf.write<size_type>((size_type)GP_map_data.size());
    for (int ti=0; ti<GP_map_data.size(); ++ti) {
        osf.writeVector(*GP_map_data.at(ti));
    }
    space->write_to_file(osf);
}

Population_checkpoint::Population_checkpoint(iSimfile& isf) {
    seed = isf.read<seed_type>();
    generation = isf.read<time_type>();
    cputime = isf.read<double>();
    pop_size = isf.read<int>();
    // Genetics:
    genetics = Genetics_sample::create_from_file(isf);
    size_type genelist_count = isf.read<size_type>(); // number of loci, really
    geneListsCopy.clear();
    geneListsCopy.reserve(genelist_count);
    for (int li=0; li<genelist_count; ++li) {
        geneListsCopy.push_back(GeneList(isf));
    }
    // Trait parameters (Omnigenic GP maps):
    size_type GP_map_data_count = isf.read<size_type>();
    GP_map_data.clear();
    for (int mi=0; mi<GP_map_data_count; ++mi) {
        GP_map_data.push_back(new std::vector<double>());
        isf.readVector<double>(*GP_map_data.back());
    }
    // Space:
    space = Space_sample::create_from_file(isf);
}


Microsample::Microsample(Population& pop, char option, double cputime) :
Sample_base(pop.getAge(), cputime, pop.size()){
    this->option = option;
    means.clear();
    covariances.clear();
    switch (option) {
        case 'n':
        break;
        case 'm':
        calc_means(pop);
        break;
        case 'v':
        calc_means(pop);
        calc_variances(pop);
        break;
        case 'c':
        calc_means(pop);
        calc_covariances(pop);
        break;
        default:
        std::cout << "Microsample error.'\n'";
        exit(1);
        break;
    }
}

void Microsample::write_to_file(oSimfile &osf) {
    osf.write<char>(option);
    osf.write<time_type>(generation);
    osf.write<double>(cputime);
    osf.write<size_type>(pop_size);
    osf.writeVector<traitType>(means);
    osf.writeVector<traitType>(covariances);
}

Microsample::Microsample(iSimfile& isf) {
    option = isf.read<char>();
    generation = isf.read<time_type>();
    cputime = isf.read<double>();
    pop_size = isf.read<size_type>();
    isf.readVector<traitType>(means);
    isf.readVector<traitType>(covariances);
}

void Microsample::calc_means(Population &pop) {
    for (Trait* tr : pop.getTraits()) {
        if (!tr->is_constant()) {
            for (int dim=0; dim<tr->get_dims(); ++dim) {
                means.push_back(tr->row_mean(dim));
            }
        }
    }
    // mean positions:
    for (int dim=0; dim<pop.getSpace().getDims(); ++dim) {
        means.push_back(pop.getSpace().get_mean(dim));
    }
}

void Microsample::calc_variances(Population& pop) {
    int mean_i = 0;
    for (Trait* tr : pop.getTraits()) {
        if (!tr->is_constant()) {
            for (int dim=0; dim<tr->get_dims(); ++dim) {
                covariances.push_back(tr->row_variance(dim,means.at(mean_i)));
                ++mean_i;
            }
        }
    }
    // spatial variation:
    for (int dim=0; dim<pop.getSpace().getDims(); ++dim) {
        covariances.push_back(pop.getSpace().get_variance(dim,means.at(mean_i)));
        ++mean_i;
    }
}

void Microsample::calc_covariances(Population& pop) {
    // Calculate (the upper triangle of) one major covariance matrix
    // The covariance between all trait dimensions and all spatial dimensions
    // First collect all data in one big matrix:
    int total_dims = pop.calc_total_data_dimensions();
    Matrix<traitType> bigM(pop_size,total_dims);
    int Mdim=0;
    for (Trait* tr : pop.getTraits()) {
        if (!tr->is_constant()) {
            for (int Tdim=0; Tdim<tr->get_dims(); ++Tdim) {
                for (int ind=0; ind<pop_size; ++ind) {
                    bigM(ind,Mdim) = (*tr)(Tdim,ind);
                }
                ++Mdim;
            }
        }
    }
    for (int Sdim=0; Sdim<pop.getSpace().getDims(); ++Sdim) {
        for (int ind=0; ind<pop_size; ++ind) {
            bigM(ind,Mdim) = pop.getSpace().getPosition(ind, Sdim);
        }
        ++Mdim;
    }
    // next, calculate all covariances, including variances
    covariances.reserve(total_dims*(total_dims+1)/2);
    covariances.clear();
    for (int i1=0; i1<total_dims; ++i1) {
        double m1 = means.at(i1);
        for (int i2=i1; i2<total_dims; ++i2) {
            double m2 = means.at(i2);
            traitType* x1 = &bigM(0,i1);
            traitType* x2 = &bigM(0,i2);
            double prodsum = 0;
            for (int ind=0; ind<pop_size; ++ind) {
                prodsum+= (*x1-m1)*(*x2-m2);
                ++x1;
                ++x2;
            }
            covariances.push_back( prodsum/pop_size );//- means.at(i1)*means.at(i2) );
        }
    }
}
