#include<iostream>
#include<vector>
#include<algorithm> // random_shuffle, swap
#include<cstdlib> // rand, srand
#include<ctime> // clock


using namespace std;


void readInput(vector<vector<double> > &matrix){
    double dist;
    for(unsigned i=0; i<matrix.size(); ++i){
        for(unsigned j=i+1; j<matrix.size(); ++j){
            cin >> dist >> dist >> dist;
            matrix[i][j] = matrix[j][i] = dist;
        }
    }
    for(unsigned i=0; i<matrix.size(); ++i)
        matrix[i][i] = 0.0;
}




struct Solution{
    vector<bool> genes;
    double fitness;
    
    Solution(){
        fitness = -1.0;
    }
    
    Solution(unsigned n){
        genes = vector<bool>(n,false);
        fitness = -1.0;
    }
};




void evaluateFitness(Solution &sol, const vector<vector<double> > &matrix){
    sol.fitness = 0.0;
    vector<unsigned> selected;
    
    for(unsigned i=0; i<sol.genes.size(); ++i){
        if(sol.genes[i])
            selected.push_back(i);
    }
    
    for(unsigned i=0; i<selected.size(); ++i){
        for(unsigned j=i+1; j<selected.size(); ++j)
            sol.fitness += matrix[selected[i]][selected[j]];
    }
}




void generateRandomSolution(Solution &sol, unsigned m){
    unsigned random_position;
    
    for(unsigned i=0; i<m; ++i){
        // genero entero entre 0 y sol.genes.size()
        random_position = rand() % sol.genes.size();
        // incremento random_position hasta que encuentre
        //    una posición no seleccionada
        while(sol.genes[random_position])
            random_position = (random_position+1) % sol.genes.size();
        sol.genes[random_position] = true;
    }
    
    random_shuffle(sol.genes.begin(), sol.genes.end());
}





void selection(const vector<Solution> &population, vector<unsigned> &selected_indivs){
    unsigned first_indiv, second_indiv;
    
    for(unsigned i=0; i<selected_indivs.size(); ++i){
        // elijo un primer individuo aleatoriamente
        first_indiv = rand() % population.size();
        // elijo un segundo individuo aleatoriamente
        second_indiv = rand() % population.size();
        // selecciono por torneo binario al que tiene mejor fitness
        if(population[first_indiv].fitness >= population[second_indiv].fitness)
            selected_indivs[i] = first_indiv;
        else
            selected_indivs[i] = second_indiv;
    }
}





double calculateContribution(const vector<unsigned> &selected, unsigned gen, const vector<vector<double> > &matrix){
    double contribution = 0.0;
    
    for(unsigned i=0; i<selected.size(); ++i)
        contribution += matrix[gen][selected[i]];
    
    return contribution;
}





void repair(Solution &sol, unsigned m, const vector<vector<double> > &matrix){
    vector<unsigned> selected;
    unsigned best_gen;
    double contrib, best_contrib;
    
    // guardo las posiciones de los genes con un 1
    for(unsigned i=0; i<sol.genes.size(); ++i){
        if(sol.genes[i])
            selected.push_back(i);
    }
    
    // mientras que la solución no sea factible por exceso de 1's
    while(selected.size() > m){
        // pongo a 0 el gen que tenga una mayor contribución en la solución
        best_contrib = -1.0;
        for(unsigned i=0; i<selected.size(); ++i){
            contrib = calculateContribution(selected,selected[i],matrix);
            if(contrib > best_contrib){
                best_gen = i;
                best_contrib = contrib;
            }
        }
        
        sol.genes[selected[best_gen]] = false;
        swap(selected[best_gen],selected.back());
        selected.pop_back();
    }
    // mientras que la solución no sea factible por falta de 1's
    while(selected.size() < m){
        // pongo a 1 el gen que aporte una mayor contribución a la solución
        best_contrib = -1.0;
        for(unsigned gen=0; gen<sol.genes.size(); ++gen){
            if(!sol.genes[gen]){
                contrib = calculateContribution(selected,gen,matrix);
                if(contrib > best_contrib){
                    best_gen = gen;
                    best_contrib = contrib;
                }
            }
        }
        
        sol.genes[best_gen] = true;
        selected.push_back(best_gen);
    }
}






void cross_uniform(const vector<Solution> &population, vector<Solution> &new_population, 
                   vector<unsigned> &selected_indivs, 
                   unsigned m, const vector<vector<double> > &matrix){
    // inicializo el fitness de los hijos a -1 (no calculado todavía)
    new_population[0].fitness = -1.0;
    new_population[1].fitness = -1.0;
    
    for(unsigned j=0; j<matrix.size(); ++j){
        // si los padres comparten el mismo gen, ambos hijos reciben ese gen
        if(population[selected_indivs[0]].genes[j] == population[selected_indivs[1]].genes[j]){
            new_population[0].genes[j] = population[selected_indivs[0]].genes[j];
            new_population[1].genes[j] = population[selected_indivs[0]].genes[j];
        } // si no, cada hijo recibe el gen de algún padre aleatoriamente
        else{
            new_population[0].genes[j] = population[selected_indivs[rand()%2]].genes[j];
            new_population[1].genes[j] = population[selected_indivs[rand()%2]].genes[j];
        }
    }
    
    // reparo a ambos hijos
    repair(new_population[0], m, matrix);
    repair(new_population[1], m, matrix);
}





void mutation(vector<Solution> &new_population, unsigned m, const vector<vector<double> > &matrix){
    for(Solution &sol : new_population){
        // muto con probabilidad 0.1
        if((rand() % 10) == 0){
            unsigned random_gene1, random_gene0, gene1_to_mutate, gene0_to_mutate;
            
            // elijo aleatoriamente dos genes de ese individuo con valores contrarios
            random_gene1 = (rand() % m) + 1;
            random_gene0 = (rand() % (matrix.size()-m)) + 1;
            unsigned i=0;
            while(random_gene1+random_gene0>0){
                if(sol.genes[i]){
                    if(random_gene1>0){
                        --random_gene1;
                        if(random_gene1==0)
                            gene1_to_mutate = i;
                    }
                }
                else{
                    if(random_gene0>0){
                        --random_gene0;
                        if(random_gene0==0)
                            gene0_to_mutate = i;
                    }
                }
                ++i;
            }
            // muto los genes elegidos
            sol.genes[gene1_to_mutate] = false;
            sol.genes[gene0_to_mutate] = true;
        }
    }
}





void replacement(vector<Solution> &population, vector<Solution> &new_population, const vector<vector<double> > &matrix){
    unsigned worst_solution;
    
    // calculo el fitness de ambos hijos
    for(Solution &sol : new_population)
        evaluateFitness(sol,matrix);
    
    // Tomo como primer hijo a aquel con mejor fitness de los dos
    if(new_population[1].fitness > new_population[0].fitness)
        swap(new_population[0],new_population[1]);
    
    // hallo la peor solución de la población actual
    worst_solution = 0;
    for(unsigned i=1; i<population.size(); ++i){
        if(population[i].fitness < population[worst_solution].fitness)
            worst_solution = i;
    }
    
    // si el primer hijo es mejor que la peor solución en la población actual, la sustituye
    if(new_population[0].fitness > population[worst_solution].fitness){
        swap(new_population[0],population[worst_solution]);
        
        // hallo la peor solución de la población actual
        worst_solution = 0;
        for(unsigned i=1; i<population.size(); ++i){
            if(population[i].fitness < population[worst_solution].fitness)
                worst_solution = i;
        }
        
        // si el segundo hijo es mejor que la peor solución en la población actual, la sustituye
        if(new_population[1].fitness > population[worst_solution].fitness)
            swap(new_population[1],population[worst_solution]);
    }
}






void AGEuniform(const vector<vector<double> > &matrix, unsigned m, unsigned population_size = 50){
    clock_t start_time, elapsed;
    start_time = clock();
    
    unsigned best_solution;
    vector<Solution> population(population_size, Solution(matrix.size()));
    vector<Solution> new_population(2, Solution(matrix.size()));
    vector<unsigned> selected_indivs(2, 0);
    unsigned evaluations;
    const unsigned MAX_EVALUATIONS = 100000;
    
    // genero población aleatoria
    for(Solution &sol : population){
        generateRandomSolution(sol, m);
        evaluateFitness(sol, matrix);
    }
    evaluations = population_size;
    
    while(evaluations < MAX_EVALUATIONS){
        // selección mediante torneo binario
        selection(population, selected_indivs);
        // cruce uniforme
        cross_uniform(population, new_population, selected_indivs, m, matrix);
        // mutación
        mutation(new_population, m, matrix);
        // reemplazamiento
        replacement(population, new_population, matrix);
        evaluations += 2;
    }
    
    // hallo la posición de la mejor solución
    best_solution = 0;
    for(unsigned i=1; i<population_size; ++i){
        if(population[i].fitness > population[best_solution].fitness)
            best_solution = i;
    }
    
    elapsed = clock() - start_time;
    
    cout << fixed;
    cout << population[best_solution].fitness << "\t" << (double) elapsed / CLOCKS_PER_SEC << endl;
}




int main(int argc, char **argv){
    unsigned n, m;
    
    cin >> n >> m;
    vector<double> v(n, 0);
    vector<vector<double> > matrix(n, v);
    readInput(matrix);
    
    srand(1);
    AGEuniform(matrix, m);
    
    return 0;
}
