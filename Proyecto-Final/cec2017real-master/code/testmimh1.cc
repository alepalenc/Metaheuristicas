extern "C" {
#include "cec17.h"
}
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm> //random_shuffle
#include <set>
#include <random>

using namespace std;




void clip(vector<double> &sol, double lower, double upper) {
  for (auto &val : sol) {
    if (val < lower) {
      val = lower;
    }
    else if (val > upper) {
      val = upper;
    }
  }
}

void increm_bias(vector<double> &bias, vector<double> dif) {
  for (unsigned i = 0; i < bias.size(); i++) {
    bias[i] = 0.2*bias[i]+0.4*(dif[i]+bias[i]);
  }
}

void decrement_bias(vector<double> &bias, vector<double> dif) {
  for (unsigned i = 0; i < bias.size(); i++) {
    bias[i] = bias[i]-0.4*(dif[i]+bias[i]);
  }
}

/**
 * Aplica el Solis Wets
 *
 * @param  sol solucion a mejorar.
 * @param fitness fitness de la solución.
 */
template <class Random>
void soliswets(vector<double> &sol, double &fitness, double delta, int maxevals, double lower, double upper, Random &random) {
    const size_t dim = sol.size();
    vector<double> bias (dim), dif (dim), newsol (dim);
    double newfit;
    size_t i;

    int evals = 0;
    int num_success = 0;
    int num_failed = 0;

    while (evals < maxevals) {
        std::uniform_real_distribution<double> distribution(0.0, delta);

        for (i = 0; i < dim; i++) {
            dif[i] = distribution(random);
            newsol[i] = sol[i] + dif[i] + bias[i];
        }

        clip(newsol, lower, upper);
        newfit = cec17_fitness(&newsol[0]);
        evals += 1;

        if (newfit < fitness) {
            sol = newsol;
            fitness = newfit;
            increm_bias(bias, dif);
            num_success += 1;
            num_failed = 0;
        }
        else if (evals < maxevals) {

            for (i = 0; i < dim; i++)
                newsol[i] = sol[i] - dif[i] - bias[i];
            
            clip(newsol, lower, upper);
            newfit = cec17_fitness(&newsol[0]);
            evals += 1;
            
            if (newfit < fitness) {
                sol = newsol;
                fitness = newfit;
                decrement_bias(bias, dif);
                num_success += 1;
                num_failed = 0;
            }
            else {
                for (i = 0; i < dim; i++)
                  bias[i] /= 2;
                
                num_success = 0;
                num_failed += 1;
            }
        }

        if (num_success >= 5) {
            num_success = 0;
            delta *= 2;
        }
        else if (num_failed >= 3) {
            num_failed = 0;
            delta /= 2;
        }
    }

}





struct Solution{
    vector<double> v;
    double fitness = -1;
    
    Solution(int dim) : v(dim){}
};



const bool VERBOSE=false;



double distance(const Solution &sol1, const Solution &sol2){
    double d = 0.0;
    for(unsigned i=0; i<sol1.v.size(); ++i)
        d += (sol1.v[i]-sol2.v[i])*(sol1.v[i]-sol2.v[i]);
    return d;
}

void findBestSolutions(vector<Solution> &population, int best_sols_to_find){
    for(unsigned i=0; i<best_sols_to_find; ++i){
        unsigned k = i;
        for(unsigned j=i+1; j<population.size(); ++j){
            if(population[k].fitness > population[j].fitness)
                k = j;
        }
        swap(population[i], population[k]);
    }
}

template <class Random>
void colonization(vector<Solution> &population, int num_childs, double lower, double upper, Random &random){
    int dim = population[0].v.size();
    std::uniform_real_distribution<double> distrib01(0.1, 1.0);
    
    random_shuffle(population.begin(), population.end());
    
    // Generación de hijos en el centro del espacio de búsqueda
    for(unsigned i=0; i<num_childs/2; ++i){
        Solution child(dim);
        for(unsigned j=0; j<dim; ++j)
            child.v[j] = (population[2*i].v[j] + population[2*i+1].v[j])/2;
        child.fitness = cec17_fitness(&child.v[0]);
        population.push_back(child);
    }
    
    // Generación de hijos en los límites del espacio de búsqueda
    for(unsigned i=num_childs; i<(3*num_childs)/2; ++i){
        Solution child(dim);
        double t = distrib01(random);
        for(unsigned j=0; j<dim; ++j){
            if(population[i].v[j] >= 0.0)
                child.v[j] = population[i].v[j] + t*(upper-population[i].v[j]);
            else
                child.v[j] = population[i].v[j] + t*(lower-population[i].v[j]);
        }
        child.fitness = cec17_fitness(&child.v[0]);
        population.push_back(child);
    }
}

template <class Random>
void development(vector<Solution> &population, int evals_per_sol, int sols_to_develop, double lower, double upper, Random &random){
    if(sols_to_develop > population.size())
        sols_to_develop = population.size();
    
    findBestSolutions(population, sols_to_develop);
    
    for(unsigned i=0; i<sols_to_develop; ++i)
        soliswets(population[i].v, population[i].fitness, 0.2, evals_per_sol, lower, upper, random);
}

void removeSolutions(vector<Solution> &population, int s, double threshold){
    set<unsigned> sols_to_remove;
    
    for(unsigned i=0; i<population.size() and sols_to_remove.size()<population.size()-s; ++i){
        for(unsigned j=i+1; j<population.size() and sols_to_remove.size()<population.size()-s; ++j){
            if(distance(population[i], population[j]) < threshold){
                if(population[i].fitness < population[j].fitness)
                    sols_to_remove.insert(j);
                else if(population[i].fitness > population[j].fitness)
                    sols_to_remove.insert(i);
            }
        }
    }
    
    if(VERBOSE)
        cout << "\tsols_to_remove.size() = " << sols_to_remove.size() << endl;
    
    set<unsigned>::reverse_iterator rit;
    for(rit=sols_to_remove.rbegin(); rit!=sols_to_remove.rend(); ++rit){
        swap(population[*rit], population.back());
        population.pop_back();
    }
}

void invasion(vector<Solution> &population, int s){
    double inv_num_distances = 2.0/(population.size()*(population.size()-1));
    double mean_distance = 0.0;
    
    // Distancia media
    for(unsigned i=0; i<population.size(); ++i){
        for(unsigned j=i+1; j<population.size(); ++j){
            mean_distance += distance(population[i], population[j])*inv_num_distances;
        }
    }
    
    random_shuffle(population.begin(), population.end());
    
    if(VERBOSE)
        cout << "\tmean_distance = " << mean_distance << endl;
    
    removeSolutions(population, s, mean_distance/16.0);
    removeSolutions(population, s, mean_distance/8.0);
    removeSolutions(population, s, mean_distance/4.0);
    removeSolutions(population, s, mean_distance/2.0);
    removeSolutions(population, s, 10000000.0);
}


/**
 * Aplica mi metaheurística
 *
 * @param fitness fitness de la solución final.
 */
template <class Random>
void mimh1(double &fitness, int max_evals, int dim, double lower, double upper, Random &random){
    vector<Solution> population;
    Solution random_sol(dim);
    std::uniform_real_distribution<double> distribution(lower, upper);
    
    int s_init = 1750+25*dim;//5*(90-dim)*dim/2;
    int s = s_init;
    int fraction = 9+dim/10;//11-dim/10;
    int p;
    int sols_to_develop;
    int development_evals = 500;
    int evals = 0;
    
    for(unsigned i=0; i<s_init; ++i){
        for(unsigned j=0; j<dim; ++j)
            random_sol.v[j] = distribution(random);
        random_sol.fitness = cec17_fitness(&random_sol.v[0]);
        population.push_back(random_sol);
    }
    evals += s_init;
    if(VERBOSE){
        cout << "iteration: " << endl;
        cout << "\ts = " << s << endl;
        cout << "\tevals = " << evals << endl;
        cout << "\tpopulation.size() = " << population.size() << endl;
        int k_aux = 0;
        for(unsigned i=0; i<population.size(); ++i){
            if(population[i].fitness < population[k_aux].fitness)
                k_aux = i;
        }
        cout << "\tbest_sol.fitness = " << population[k_aux].fitness << endl;
    }
    
    while(evals < max_evals-s_init and s > 20){
        p = ceil(s/fraction);
        p += p%2; // para que p sea par
        colonization(population, p, lower, upper, random);
        evals += p;
        sols_to_develop = floor((s_init-p)/development_evals);
        development(population, development_evals, sols_to_develop, lower, upper, random);
        if(VERBOSE){
            cout << "iteration (COL-DEV-INV): " << endl;
            cout << "\ts = " << s << endl;
            cout << "\tp = " << p << endl;
            cout << "\tevals_ls = " << development_evals*sols_to_develop << endl;
        }
        evals += development_evals*sols_to_develop;
        s -= p;
        invasion(population, s);
        if(VERBOSE){
            cout << "\tevals = " << evals << endl;
            cout << "\tpopulation.size() = " << population.size() << endl;
            int k_aux = 0;
            for(unsigned i_aux=0; i_aux<population.size(); ++i_aux){
                if(population[i_aux].fitness < population[k_aux].fitness)
                    k_aux = i_aux;
            }
            cout << "\tbest_sol.fitness = " << population[k_aux].fitness << endl;
        }
    }
    
    sols_to_develop = 5;
    findBestSolutions(population, sols_to_develop);
    
    while(evals < max_evals-s_init){
        for(unsigned i=0; i<sols_to_develop; ++i)
            soliswets(population[i].v, population[i].fitness, 0.2, development_evals, lower, upper, random);
        evals += development_evals*sols_to_develop;
        if(VERBOSE){
            cout << "iteration: " << endl;
            cout << "\tevals_ls = " << development_evals*sols_to_develop << endl;
            cout << "\tevals = " << evals << endl;
            cout << "\tpopulation.size() = " << population.size() << endl;
            int k_aux = 0;
            for(unsigned i_aux=0; i_aux<population.size(); ++i_aux){
                if(population[i_aux].fitness < population[k_aux].fitness)
                    k_aux = i_aux;
            }
            cout << "\tbest_sol.fitness = " << population[k_aux].fitness << endl;
        }
    }
    
    findBestSolutions(population, 1);
    soliswets(population[0].v, population[0].fitness, 0.2, max_evals-evals, lower, upper, random);
    
    fitness = population[0].fitness;
}



int main(){
    int dim = 10;
    int range_seed = 10;
    double upper = 100.0;
    double lower = -100.0;
    std::uniform_real_distribution<> dis(lower, upper);
    
    for(int seed=1; seed<=range_seed; ++seed){
        cout << "///// SEED " << seed << " /////" << endl;
        for(int funcid=1; funcid<=30; ++funcid){
            double fitness;
            
            cec17_init("mimh1", funcid, dim);
            
            cerr <<"Warning: output by console, if you want to create the output file you have to comment cec17_print_output()" <<endl;
            cec17_print_output(); // Comment to generate the output file
            
            std::mt19937 gen(seed); // Inicio semilla
            
            mimh1(fitness, 10000*dim, dim, lower, upper, gen);
            
            double error = cec17_error(fitness);
            cout << "Error[F" << funcid << "]: " << scientific << error << endl;
        }
    }
    
    return 0;
}
