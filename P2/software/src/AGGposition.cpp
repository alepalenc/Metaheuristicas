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





void cross_position(vector<Solution> &population, vector<Solution> &new_population, 
                    vector<unsigned> &selected_indivs, 
                    const vector<vector<double> > &matrix){
    unsigned couples_to_cross = (unsigned)((population.size()/2) * 0.7); // número de parejas a cruzar
    vector<bool> first_father_remains;
    vector<unsigned> positions_to_complete;
    
    // cruzo parejas
    for(unsigned i=0; i<couples_to_cross; ++i){
        // inicializo el fitness de los hijos a -1 (no calculado todavía)
        new_population[2*i].fitness = -1.0;
        new_population[2*i+1].fitness = -1.0;
        
        // vacío los vectores de restos del primer padre y de posiciones a completar
        first_father_remains.clear();
        positions_to_complete.clear();
        for(unsigned j=0; j<matrix.size(); ++j){
            // si los padres comparten el mismo gen, ambos hijos reciben ese gen
            if(population[selected_indivs[2*i]].genes[j] == population[selected_indivs[2*i+1]].genes[j]){
                new_population[2*i].genes[j] = population[selected_indivs[2*i]].genes[j];
                new_population[2*i+1].genes[j] = population[selected_indivs[2*i]].genes[j];
            } // si no, guardo el gen del primer padre y la posición en la que se encuentra
            else{
                first_father_remains.push_back(population[selected_indivs[2*i]].genes[j]);
                positions_to_complete.push_back(j);
            }
        }
        // barajo el vector de genes restantes del primer padre
        random_shuffle(first_father_remains.begin(), first_father_remains.end());
        
        // completo el primer hijo con los genes restantes del primer padre barajados
        for(unsigned j=0; j<first_father_remains.size(); ++j)
            new_population[2*i].genes[positions_to_complete[j]] = first_father_remains[j];
        
        // vuelvo a barajar el vector de genes restantes del primer padre
        random_shuffle(first_father_remains.begin(), first_father_remains.end());
        
        // completo el segundo hijo con los genes restantes del primer padre barajados
        for(unsigned j=0; j<first_father_remains.size(); ++j)
            new_population[2*i+1].genes[positions_to_complete[j]] = first_father_remains[j];
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






void AGGposition(const vector<vector<double> > &matrix, unsigned m, unsigned population_size = 50){
    clock_t start_time, elapsed;
    start_time = clock();
    
    unsigned best_solution;
    vector<Solution> population(population_size, Solution(matrix.size()));
    vector<Solution> new_population(population_size, Solution(matrix.size()));
    vector<unsigned> selected_indivs(population_size, 0);
    unsigned evaluations;
    const unsigned MAX_EVALUATIONS = 100000;
    
    // genero población aleatoria
    for(Solution &sol : population){
        generateRandomSolution(sol, m);
        evaluateFitness(sol, matrix);
    }
    evaluations = population_size;
    
    // hallo la posición de la mejor solución
    best_solution = 0;
    for(unsigned i=0; i<population_size; ++i){
        if(population[i].fitness > population[best_solution].fitness)
            best_solution = i;
    }
    
    while(evaluations < MAX_EVALUATIONS){
        // selección mediante torneo binario
        selection(population, selected_indivs);
        // cruce posición
        cross_position(population, new_population, selected_indivs, matrix);
        // mutación
        mutation(new_population, m, matrix);
        // reemplazamiento
        evaluations += replacement(population, new_population, best_solution, matrix);
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
    AGGposition(matrix, m);
    
    return 0;
}
