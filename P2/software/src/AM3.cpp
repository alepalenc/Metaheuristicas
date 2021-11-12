#include<iostream>
#include<vector>
#include<algorithm> // swap
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





double calculateContribution(const vector<unsigned> &selected, unsigned gen, const vector<vector<double> > &matrix){
    double contribution = 0.0;
    
    for(unsigned i=0; i<selected.size(); ++i)
        contribution += matrix[gen][selected[i]];
    
    return contribution;
}






unsigned localSearch(Solution &sol, const vector<vector<double> > &matrix){
    vector<unsigned> selected, not_selected;
    double element_contrib, worst_element_contrib;
    bool exchange;
    unsigned evaluations;
    const unsigned MAX_EVALUATIONS = 400;
    
    // pasar la solución en codificación binaria a una codificación en conjuntos de enteros
    //   selected: conjunto de elementos seleccionados
    //   not_selected: conjunto de elementos no seleccionados
    for(unsigned i=0; i<sol.genes.size(); ++i){
        if(sol.genes[i])
            selected.push_back(i);
        else
            not_selected.push_back(i);
    }
    
    evaluations = 0;
    do{
        exchange = false;
        // fijo un elemento (posición de un gen) en la solución
        for(unsigned i=0; i<selected.size() and !exchange and evaluations<MAX_EVALUATIONS; ++i){
            // obtengo el siguiente elemento seleccionado que tiene una menor contribución
            //   y lo coloco en la posición i-ésima
            worst_element_contrib = calculateContribution(selected, selected[i], matrix);
            for(unsigned j=i+1; j<selected.size(); ++j){
                element_contrib = calculateContribution(selected, selected[j], matrix);
                if(element_contrib < worst_element_contrib){
                    swap(selected[i], selected[j]);
                    worst_element_contrib = element_contrib;
                }
            }
            
            // fijo una posición del conjunto de genes no seleccionados
            for(unsigned j=0; j<not_selected.size() and !exchange and evaluations<MAX_EVALUATIONS; ++j){
                element_contrib = calculateContribution(selected, not_selected[j], matrix) - matrix[selected[i]][not_selected[j]];
                ++evaluations;
                // si la contribución del elemento no seleccionado es mejor que la del seleccionado, los intercambio
                if(element_contrib > worst_element_contrib){
                    swap(selected[i], not_selected[j]);
                    exchange = true;
                }
            }
        }
    }while(exchange && evaluations < MAX_EVALUATIONS);
    
    // actualizo los genes de la solución en codificación binaria
    for(unsigned i=0; i<selected.size(); ++i)
        sol.genes[selected[i]] = true;
    for(unsigned i=0; i<not_selected.size(); ++i)
        sol.genes[not_selected[i]] = false;
        
    // actualizo el fitness de la solución
    sol.fitness = 0.0;
    for(unsigned i=0; i<selected.size(); ++i){
        for(unsigned j=i+1; j<selected.size(); ++j)
            sol.fitness += matrix[selected[i]][selected[j]];
    }
    
    return evaluations;
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
    unsigned couples_to_cross = (unsigned)((population.size()/2) * 0.7); // número de parejas a cruzar
    
    // cruzo parejas
    for(unsigned i=0; i<couples_to_cross; ++i){
        // inicializo el fitness de los hijos a -1 (no calculado todavía)
        new_population[2*i].fitness = -1.0;
        new_population[2*i+1].fitness = -1.0;
        
        for(unsigned j=0; j<matrix.size(); ++j){
            // si los padres comparten el mismo gen, ambos hijos reciben ese gen
            if(population[selected_indivs[2*i]].genes[j] == population[selected_indivs[2*i+1]].genes[j]){
                new_population[2*i].genes[j] = population[selected_indivs[2*i]].genes[j];
                new_population[2*i+1].genes[j] = population[selected_indivs[2*i]].genes[j];
            } // si no, cada hijo recibe el gen de algún padre aleatoriamente
            else{
                new_population[2*i].genes[j] = population[selected_indivs[2*i+(rand()%2)]].genes[j];
                new_population[2*i+1].genes[j] = population[selected_indivs[2*i+(rand()%2)]].genes[j];
            }
        }
        
        // reparo a ambos hijos
        repair(new_population[2*i], m, matrix);
        repair(new_population[2*i+1], m, matrix);
    }
    
    // el resto de individuos pasar a la nueva población sin cruzar
    for(unsigned i=2*couples_to_cross; i<selected_indivs.size(); ++i)
        new_population[i] = population[selected_indivs[i]];
}





void mutation(vector<Solution> &new_population, unsigned m, const vector<vector<double> > &matrix){
    unsigned genes_to_mutate = (unsigned)(new_population.size() * 0.1); // número de genes a mutar
    unsigned random_indiv, random_gene1, random_gene0, gene1_to_mutate, gene0_to_mutate;
    
    for(unsigned k=0; k<genes_to_mutate; ++k){
        // elijo aleatoriamente un individuo
        random_indiv = rand() % new_population.size();
        // elijo aleatoriamente dos genes de ese individuo con valores contrarios
        random_gene1 = (rand() % m) + 1;
        random_gene0 = (rand() % (matrix.size()-m)) + 1;
        unsigned i=0;
        while(random_gene1+random_gene0>0){
            if(new_population[random_indiv].genes[i]){
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
        new_population[random_indiv].genes[gene1_to_mutate] = false;
        new_population[random_indiv].genes[gene0_to_mutate] = true;
    }
}





unsigned replacement(vector<Solution> &population, vector<Solution> &new_population, 
                     unsigned &best_solution_population, const vector<vector<double> > &matrix){
    unsigned best_solution_new_population, worst_solution_new_population, inc_evaluations;
    
    // calculo el fitness de todas las soluciones de la nueva población que no lo tengan calculado
    inc_evaluations = 0;
    for(Solution &sol : new_population){
        if(sol.fitness < 0.0){
            evaluateFitness(sol,matrix);
            ++inc_evaluations;
        }
    }
    
    // hallo la posición del mejor individuo de la nueva población
    best_solution_new_population = 0;
    for(unsigned i=1; i<new_population.size(); ++i){
        if(new_population[i].fitness > new_population[best_solution_new_population].fitness)
            best_solution_new_population = i;
    }
    
    // si la mejor solución de la anterior población es mejor que la mejor solución de la nueva población
    if(new_population[best_solution_new_population].fitness < population[best_solution_population].fitness){
        // calculo la posición de la peor solución de la nueva población
        worst_solution_new_population = 0;
        for(unsigned i=1; i<new_population.size(); ++i){
            if(new_population[i].fitness < new_population[worst_solution_new_population].fitness)
                worst_solution_new_population = i;
        }
        
        // sustituyo la peor solución de la nueva población por la mejor solución de la anterior población
        swap(new_population[worst_solution_new_population],population[best_solution_population]);
        // actualizo la variable que guarda la posición de la mejor solución de la población actual
        best_solution_population = worst_solution_new_population;
    }
    else{
        // actualizo la variable que guarda la posición de la mejor solución de la población actual
        best_solution_population = best_solution_new_population;
    }
    
    // sustituyo a la anterior población por la nueva población
    swap(population,new_population);
    
    return inc_evaluations;
}






void AM3(const vector<vector<double> > &matrix, unsigned m, unsigned population_size = 50){
    clock_t start_time, elapsed;
    start_time = clock();
    
    unsigned best_solution;
    vector<Solution> population(population_size, Solution(matrix.size()));
    vector<Solution> new_population(population_size, Solution(matrix.size()));
    vector<unsigned> selected_indivs(population_size, 0);
    unsigned evaluations;
    const unsigned MAX_EVALUATIONS = 100000, GENERATIONS_BETWEEN_LOCAL_SEARCH = 10;
    const unsigned SOLUTIONS_TO_APPLY_LOCAL_SEARCH = population_size*0.1;
    
    // genero población aleatoria
    for(Solution &sol : population){
        generateRandomSolution(sol, m);
        evaluateFitness(sol, matrix);
    }
    evaluations = population_size;
    
    // hallo la posición de la mejor solución
    best_solution = 0;
    for(unsigned i=1; i<population_size; ++i){
        if(population[i].fitness > population[best_solution].fitness)
            best_solution = i;
    }
    
    while(evaluations < MAX_EVALUATIONS){
        for(unsigned k=0; k<GENERATIONS_BETWEEN_LOCAL_SEARCH and evaluations < MAX_EVALUATIONS; ++k){
            // selección mediante torneo binario
            selection(population, selected_indivs);
            // cruce uniforme con reparación
            cross_uniform(population, new_population, selected_indivs, m, matrix);
            // mutación
            mutation(new_population, m, matrix);
            // reemplazamiento
            evaluations += replacement(population, new_population, best_solution, matrix);
        }
        
        // coloco las mejores soluciones al principio del vector population
        swap(population[best_solution], population[0]);
        for(unsigned i=1; i<SOLUTIONS_TO_APPLY_LOCAL_SEARCH; ++i){
            best_solution = i;
            for(unsigned j=i+1; j<population.size(); ++j){
                if(population[best_solution].fitness < population[j].fitness)
                    best_solution = j;
            }
            swap(population[i], population[best_solution]);
        }
        best_solution = 0;
        
        // aplico búsqueda local a las 10 mejores soluciones
        for(unsigned i=0; i<SOLUTIONS_TO_APPLY_LOCAL_SEARCH and evaluations<MAX_EVALUATIONS; ++i){
            evaluations += localSearch(population[i], matrix);
            if(population[i].fitness > population[best_solution].fitness)
                best_solution = i;
        }
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
    AM3(matrix, m);
    
    return 0;
}
