#include<iostream>
#include<vector>
#include<algorithm> // swap, random_shuffle
#include<cstdlib> // rand, srand
#include<ctime> // clock
#include<cassert> // assert


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
    vector<unsigned> selected;
    vector<unsigned> not_selected;
    double fitness = -1.0;
};




void evaluateFitness(Solution &sol, const vector<vector<double> > &matrix){
    sol.fitness = 0.0;
    
    for(unsigned i=0; i<sol.selected.size(); ++i){
        for(unsigned j=i+1; j<sol.selected.size(); ++j)
            sol.fitness += matrix[sol.selected[i]][sol.selected[j]];
    }
}




void generateRandomSolution(Solution &sol, unsigned n, unsigned m){
    unsigned random_element;
    
    sol = Solution();
    
    for(unsigned i=0; i<n; ++i)
        sol.not_selected.push_back(i);
    
    random_shuffle(sol.not_selected.begin(), sol.not_selected.end());
    
    for(unsigned i=0; i<m; ++i){
        // genero entero entre 0 y sol.not_selected.size()-1
        random_element = rand() % sol.not_selected.size();
        // intercambio los valores de sol.not_selected.back() y sol.not_selected[random_element]
        swap(sol.not_selected.back(), sol.not_selected[random_element]);
        // añado como gen el último elemento de not_selected
        sol.selected.push_back(sol.not_selected.back());
        // elimino el último elemento de not_selected
        sol.not_selected.pop_back();
    }
}





double contribution(const Solution &sol, unsigned element, const vector<vector<double> > &matrix){
    double contrib = 0.0;
    
    for(unsigned i=0; i<sol.selected.size(); ++i)
        contrib += matrix[element][sol.selected[i]];
    
    return contrib;
}






void localSearch(Solution &sol, unsigned max_evaluations, const vector<vector<double> > &matrix){
    double element_contrib, worst_element_contrib;
    bool exchange;
    unsigned evaluations;
    
    evaluations = 0;
    do{
        exchange = false;
        // fijo un elemento (posición de un gen) en la solución
        for(unsigned i=0; i<sol.selected.size() and !exchange and evaluations<max_evaluations; ++i){
            // obtengo el siguiente elemento seleccionado que tiene una menor contribución
            //   y lo coloco en la posición i-ésima
            worst_element_contrib = contribution(sol, sol.selected[i], matrix);
            for(unsigned j=i+1; j<sol.selected.size(); ++j){
                element_contrib = contribution(sol, sol.selected[j], matrix);
                if(element_contrib < worst_element_contrib){
                    swap(sol.selected[i], sol.selected[j]);
                    worst_element_contrib = element_contrib;
                }
            }
            
            // fijo una posición del conjunto de genes no seleccionados
            for(unsigned j=0; j<sol.not_selected.size() and !exchange and evaluations<max_evaluations; ++j){
                element_contrib = contribution(sol, sol.not_selected[j], matrix) - matrix[sol.selected[i]][sol.not_selected[j]];
                ++evaluations;
                // si la contribución del elemento no seleccionado es mejor que la del seleccionado, los intercambio
                if(element_contrib > worst_element_contrib){
                    swap(sol.selected[i], sol.not_selected[j]);
                    sol.fitness += element_contrib - worst_element_contrib;
                    exchange = true;
                }
            }
        }
    }while(exchange && evaluations < max_evaluations);
}






void mutation(Solution &sol, unsigned elements_to_mutate, const vector<vector<double> > &matrix){
    random_shuffle(sol.selected.begin(), sol.selected.end());
    random_shuffle(sol.not_selected.begin(), sol.not_selected.end());
    
    for(unsigned i=0; i<elements_to_mutate; ++i){
        sol.fitness += ( contribution(sol,sol.not_selected[i],matrix)
                       - matrix[sol.selected[i]][sol.not_selected[i]] )
                       - contribution(sol,sol.selected[i],matrix);
        swap(sol.selected[i], sol.not_selected[i]);
    }
}





void ILS(const vector<vector<double> > &matrix, unsigned m){
    clock_t start_time, elapsed;
    start_time = clock();
    
    const unsigned ITERATIONS = 10,
                   MAX_EVALUATIONS_LS = 10000,
                   ELEMENTS_TO_MUTATE = m/10;
    
    Solution sol, best_sol;
    
    generateRandomSolution(sol, matrix.size(), m);
    evaluateFitness(sol, matrix);
    localSearch(sol, MAX_EVALUATIONS_LS, matrix);
    best_sol = sol;
    
    for(unsigned i=1; i<ITERATIONS; ++i){
        mutation(sol, ELEMENTS_TO_MUTATE, matrix);
        localSearch(sol, MAX_EVALUATIONS_LS, matrix);
        if( sol.fitness > best_sol.fitness )
            best_sol = sol;
        else
            sol = best_sol;
    }
    
    elapsed = clock() - start_time;
    
    cout << fixed;
    cout << best_sol.fitness << "\t" << (double) elapsed / CLOCKS_PER_SEC << endl;
}




int main(int argc, char **argv){
    unsigned n, m;
    
    cin >> n >> m;
    vector<double> v(n, 0);
    vector<vector<double> > matrix(n, v);
    readInput(matrix);
    
    srand(1);
    ILS(matrix, m);
    
    return 0;
}
