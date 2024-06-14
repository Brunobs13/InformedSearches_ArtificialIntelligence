#Bruno Ricardo de Sá Ferreire 
#Engenharia Informática
#Nestudante 2201529

import heapq #para manipular filas de prioridade, heaps, usefull for bfs and Astar 
import itertools #iteradores que geram combinacoes de movimentos
import random#geracao de numeros aleatorios
import time #controle e manipulacao de tempo
import math #arredondamento e possiveis calculos matematicos
from mapas import instancias, coordenadas_especificas_lista #matrizes propostas
#from Algortimos import hill_climbing, best_first_search, genetic_algorithm

#distancia chebyshev, distancia mais curta com diagonais incluidas. 
def chebyshev_distance(x1, y1, x2, y2):#abs para numeros absolutos max, para o maior entre os dois 
	return max(abs(x1 - x2), abs(y1 - y2))# retornamos o maior numero que resulta das diferenças

def travel_cost(distance):#calcular o custo da viagem com base da distancia fornecida
	if distance == 0:
		return 0
	elif distance == 1:
		return 0
	elif distance == 2:
		return 1
	elif distance == 3:
		return 2
	elif distance == 4:
		return 4
	elif distance == 5:
		return 8
	else:
		return 10
	
def dentro_do_mapa(x, y, matriz):#verificar se as coordenadas estao dentro do mapa 
	return 0 <= x < len(matriz) and 0 <= y < len(matriz[0]) #x maior ou igual a prmeira coordenada e menor que as dimensões da matriz, o mesmo para o y

def gerar_posicoes_aleatorias(matriz, num_estacoes):#geramos posicoes para a matriz, passamos matriz para sabermos as dimensoes e posicoes a gerar, e num de estacoes
	linhas = len(matriz)#aqui passamos o numero de linhas para linhas
	colunas = len(matriz[0]) if linhas > 0 else 0 #se houver pelo menos uma linha ,definimos o numero de colunas
	posicoes = set()#armazenar elementos unicos com set, colecao nao ordenada,nao permite duplicatas, permite inserir e remover. Nao sera possivel gerar duas (x,y) iguais
	while len(posicoes) < num_estacoes:#loop, enquanto num_estacoes pretendido for maior que posicoes a serem atrivuidas 
		x = random.randint(0, linhas - 1)#geração aleatoria para as linhas, a partir de zero ate a ultima posicao, linhas-1 = numero total de linhas
		y = random.randint(0, colunas - 1)#o mesmo para y
		posicoes.add((x, y))  #adiciona as coordenadas ao conjunto set de posicoes para evitar repetições 
	return list(posicoes)#converte posicoes numa lista, lista pode ser usada para fazer outras conversoes, boa pratica 

def usar_coordenadas_especificas(idx):
	if idx < len(coordenadas_especificas_lista):
		return coordenadas_especificas_lista[idx]
	else:
		return []


def heuristic(map, stations):#funcao heuristica com dois argumentos, matriz com as populacoes e lista com coordenadas
	total_penalty = 0#variavel de penalizacao tendo em conta distancia e custo
	total_coverage_benefit = 0
	for i in range(len(map)):#iterar sobre a matriz com i,j
		for j in range(len(map[0])):
			if map[i][j] > 0:#verifica se posicao atual tem populacao
				min_station_dist = min(chebyshev_distance(i, j, sx, sy) for sx, sy in stations)#usa funcao cheby e calcula distancia minima entre posicao atual e estacao
				if min_station_dist <= 1:#se a distancia e 0 ou 1 
				    total_coverage_benefit += map[i][j] #adicionamos o beneficio de termos familias sem custo de deslocacao
				total_penalty += map[i][j] * travel_cost(min_station_dist)#penalizamos o custo de deslocacao
	return max(0, total_penalty - total_coverage_benefit*0.3)#retornamos a penalidade total, menos as familias no raio 1, com um peso de 30% 


def calculate_travel_cost(map, stations):#vamos calcular o custo médio de deslocamento
	total_cost = 0#custo total do deslocamento
	total_population = 0#populacao total
	for i in range(len(map)):#iterar mapa
		for j in range(len(map[0])):
			if map[i][j] > 0:#verifica se ha populacao
				closest_station_distance = min(chebyshev_distance(i, j, sx, sy) for sx, sy in stations)#calcula menor distancia entre celula atual e estaca0
				total_cost += travel_cost(closest_station_distance) * map[i][j]#total cost de deslocamento de celula atual para estacao
				total_population += map[i][j]#adiciona as familias na celula atual a variavel total_popu
	return total_cost / total_population if total_population else float('inf')#calcula custo total a dividir por numero de familias 

def minimize_cost(map, stations):
	num_stations = len(stations)
	avg_travel_cost = calculate_travel_cost(map, stations)
	return 1000 * num_stations + int(avg_travel_cost * 100)

# Algoritmo Escalada do Monte
def hill_climbing(map, num_estacoes):
	start_time = time.time()
	
	def generate_neighbors(stations):
		neighbors = []
		for i in range(len(stations)):
			x, y = stations[i]
			for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
				new_x, new_y = x + dx, y + dy
				if dentro_do_mapa(new_x, new_y, map):
					neighbor = stations[:i] + [(new_x, new_y)] + stations[i+1:]
					neighbors.append(neighbor)
		return neighbors
	
	current_solution = gerar_posicoes_aleatorias(map, num_estacoes)
	current_cost = calculate_travel_cost(map, current_solution)
	counter = {'estados_gerados': 1, 'expansoes': 0}
	
	while True:
		if (time.time() - start_time) > 60 or counter['estados_gerados'] > 100000:
				break
		neighbors = generate_neighbors(current_solution)
		neighbors_costs = [(neighbor, calculate_travel_cost(map, neighbor)) for neighbor in neighbors]
		counter['expansoes'] += len(neighbors_costs)
		neighbors_costs.sort(key=lambda x: x[1])
		
		if neighbors_costs and neighbors_costs[0][1] < current_cost:
			current_solution, current_cost = neighbors_costs[0]
			counter['estados_gerados'] += 1
		else:
			break
		
	return current_solution, current_cost, counter



# Algoritmo Melhor Primeiro
def best_first_search(map, initial_stations, limit=1000):
	start_time = time.time()
	if any(not dentro_do_mapa(x, y, map) for x, y in initial_stations):
		print("Erro: Uma ou mais estações iniciais estão fora dos limites do mapa.")
		return None, None
	
	frontier = [(heuristic(map, initial_stations), tuple(initial_stations), 0)]
	visited = set([tuple(initial_stations)])
	counter = {'estados_gerados': 1, 'expansoes': 0}
	
	
	while frontier:
		if (time.time() - start_time) > 60 or counter['estados_gerados'] > 100000:
				break
		h_cost, current_stations, depth = frontier.pop(0)
		counter['expansoes'] += 1
		
		if depth >= limit:
			continue
		
		current_travel_cost = calculate_travel_cost(map, current_stations)
		if current_travel_cost < 3:
			total_cost = minimize_cost(map, current_stations)
			stats = {
				'estados_gerados': counter['estados_gerados'],
				'expansoes': counter['expansoes'],
				'custo_total': total_cost,
				'custo_medio': math.floor(current_travel_cost * 100) / 100.0,
				'tempo_execucao': None  # This will be calculated in the calling function
			}
			return current_stations, stats
		
		successors = expand_and_generate_successors(map, current_stations)
		for successor in successors:
			if tuple(successor) not in visited:
				visited.add(tuple(successor))
				
				frontier.append((heuristic(map, successor), tuple(successor), depth + 1))
				counter['estados_gerados'] += 1
				
		frontier.sort(key=lambda x: x[0])#ordenacao logica do algortimo BFS
		
	print(f"Nenhuma solução válida encontrada após gerar {counter['estados_gerados']} estados.")
	print(f"Total de expansões: {counter['expansoes']}")
	return None, None


# Algoritmo Genético
def genetic_algorithm(map, population_size, mutation_rate, num_estacoes):
	start_time = time.time()
	
	def crossover(parent1, parent2):#combina partes de dois pais para criar um filho
		if len(parent1) > 1:
			point = random.randint(1, len(parent1) - 1)
			child = parent1[:point] + parent2[point:]
		else:
			child = parent1
		return child
	
	def mutate(solution):#muta aleatoriamente uma estação na solução
		idx = random.randint(0, len(solution) - 1)
		x, y = solution[idx]
		new_x = random.randint(0, len(map) - 1)
		new_y = random.randint(0, len(map[0]) - 1)
		solution[idx] = (new_x, new_y)
		
	def generate_initial_population():
		return [gerar_posicoes_aleatorias(map, num_estacoes) for _ in range(population_size)]
	
	population = generate_initial_population()
	generations = 0
	counter = {'estados_gerados': population_size, 'expansoes': 0}
	best_solution = None
	best_cost = float('inf')
	
	while generations < 100:
		if (time.time() - start_time) > 60 or counter['estados_gerados'] > 100000:
				break
		population_costs = [(solution, calculate_travel_cost(map, solution)) for solution in population]
		population_costs.sort(key=lambda x: x[1])
		counter['expansoes'] += len(population_costs)
		
		if population_costs[0][1] < best_cost:
			best_solution, best_cost = population_costs[0]
			
		new_population = population_costs[:2]  # Eliticismo
		while len(new_population) < population_size:
			parent1, parent2 = random.choices(population, k=2)
			child = crossover(parent1, parent2)
			if random.random() < mutation_rate:
				mutate(child)
			new_population.append((child, calculate_travel_cost(map, child)))
			counter['estados_gerados'] += 1
			
		population = [solution for solution, cost in new_population]
		generations += 1
		avg_cost =  calculate_travel_cost(map, best_solution)
		
	return best_solution, {'estados_gerados': counter['estados_gerados'], 'expansoes': counter['expansoes'], 'custo_total': minimize_cost(map, best_solution), 'custo_medio': math.floor(avg_cost * 100) / 100.0, 'tempo_execucao': None}

# Algoritmo A*
def a_star_algorithm(map, initial_stations):#recebe matriz e lista de coordenadas
	start_time = time.time()#marca o tempo de inicio de execução do algoritmo
	if any(not dentro_do_mapa(x, y, map) for x, y in initial_stations):#verifica se estacao inicial esta fora dos limites do mapa
		print("Erro: Uma ou mais estações iniciais estão fora dos limites do mapa.")
		return None, None
	initial_cost = minimize_cost(map, initial_stations)#calcular o custo total com coordenadas das estacoes iniciais 
	frontier = [(initial_cost + heuristic(map, initial_stations), tuple(initial_stations), 0, initial_cost)]#vamos ordenar os estados as serem expandidos, atraves de uma fila de prioridades, tendo em conta o custo total estimado f(n), cada item e uma tupla, expande o estado com o menor f(n)
	visited = set([tuple(initial_stations)])#vamos rastrear os estados ja visitados, set para avaliar duplicatas
	counter = {'estados_gerados': 1, 'expansoes': 0}#dicionario que armazena os estados gerados e expansoes
	while frontier:#loop principal, enquanto houver estados em frontier
		if (time.time() - start_time) > 60 or counter['estados_gerados'] > 100000:
				break#condicao de parada, tempo ou 100000 avaliacoes 
		f_cost, current_stations, g_cost, h_cost = heapq.heappop(frontier)#implementa heap em frontier e retorna menor elemento do heap, frontier e desempacatado apos aplicacao de heap 
		counter['expansoes'] += 1#aumenta o contador de expansoes
		current_travel_cost = calculate_travel_cost(map, current_stations)#calcula o custo medio de deslocamento
		if current_travel_cost < 3:#verifica se o custo medio de deslocamento e menor que 3
			total_cost = minimize_cost(map, current_stations)#aplica a funcao de custo minimizacao de custo
			stats = {#dicionario stats atualiza estisticas
				'estados_gerados': counter['estados_gerados'],
				'expansoes': counter['expansoes'],
				'custo_total': total_cost,
				'custo_medio': math.floor(current_travel_cost * 100) / 100.0,  #arredonda para baixo com duas casas decimais 
				'tempo_execucao': None  
			}
			return current_stations, stats#retorna configuracoes atuais
		expand_and_generate_states(map, current_stations, visited, frontier, counter)#chama esta funcao para expandir o estado atual e gerar novos
	print(f"Nenhuma solução válida encontrada após gerar {counter['estados_gerados']} estados.")
	print(f"Total de expansões: {counter['expansoes']}")
	return None, None

def expand_and_generate_states(map, current_stations, visited, frontier, counter):
	all_possible_moves = []#lista que armazena todos os movimentos possiveis
	for x, y in current_stations:#curret stations e a lista de coordenadas das estacoes atuais
		valid_moves = [(x + dx, y + dy) for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)] if dentro_do_mapa(x + dx, y + dy, map)]#se esta dentro do mapa atribuimos os movimentos validos, feito atraves de desempacotamento das tuplas com os movimentos possiveis em dx, dy
		valid_moves.insert(0, (x, y))#correcao apos verificar que nao gera estados a partir da primeira coordenada
		all_possible_moves.append(valid_moves if valid_moves else [(x, y)])#retorna movimentos validos, caso contrario a posicao original 
	product_of_moves = itertools.product(*all_possible_moves)#gera todas as combinacoes de movimentacoes possiveis para todas as estacoes
	for new_stations in product_of_moves:#itera sobre cada combinacao de novas estacoes
		if len(set(new_stations)) == len(new_stations): #verifica se nao ha duplicatas
			new_stations_tuple = tuple(new_stations)#converte newstation em tupla para fins de compatibilidade 
			if new_stations_tuple not in visited:# verifica se o estado new_stations_tuple ainda nao foi visitado
				visited.add(new_stations_tuple)#marca como visitado 
				new_g = minimize_cost(map, new_stations)#faz os calculos de custo
				new_h = heuristic(map, new_stations)
				new_f = new_g + new_h
				heapq.heappush(frontier, (new_f, new_stations_tuple, new_g, new_h))# adiciona o novo estado a fronteira frontier, ordenado pelo custo total f.
				counter['estados_gerados'] += 1
				
			
def expand_and_generate_successors(map, current_stations):
	all_possible_moves = []
	for x, y in current_stations:
		valid_moves = [(x + dx, y + dy) for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)] if dentro_do_mapa(x + dx, y + dy, map)]
		valid_moves.insert(0, (x, y))
		all_possible_moves.append(valid_moves if valid_moves else [(x, y)])
	return [list(new_stations) for new_stations in itertools.product(*all_possible_moves)]

#Funcao para imprimir os resultados 
def print_matrix_and_stats(matriz, stations, stats, idx):
	print(f"\nMatriz {idx + 1}:")
	station_set = set(stations)
	for i in range(len(matriz)):
		for j in range(len(matriz[0])):
			if (i, j) in station_set:
				print(f"{matriz[i][j]}#", end='  ')
			else:
				print(f"{matriz[i][j]:2}", end='  ')
		print()  #nova linha para cada matriz
		
	print("\nEstatísticas:")
	print(f"Total de estados gerados: {stats['estados_gerados']}")
	print(f"Total de expansões: {stats['expansoes']}")
	print(f"Custo total: {stats['custo_total']}")
	print(f"Custo médio: {stats['custo_medio']:.2f}")
	print(f"Tempo de execução (ms): {stats['tempo_execucao']}")
	
#processar matrizes com o algoritmo escolhido
def processar_matrizes(algoritmo, instancias, num_estacoes_por_matriz, coordenadas_especificas):
	for idx, (mapa, num_estacoes) in enumerate(zip(instancias, num_estacoes_por_matriz)):
		if num_estacoes == 0:
				print(f"Erro: Nenhuma estação inicial fornecida para a matriz {idx + 1}.")
				continue
		start_time = time.time()
		stations, stats = None, None
		
		if coordenadas_especificas:
			stations = usar_coordenadas_especificas(idx)
		else:
			stations = gerar_posicoes_aleatorias(mapa, num_estacoes)
		
		if algoritmo == 'A*':
			stations, stats = a_star_algorithm(mapa, stations)
		elif algoritmo == 'Genético':
			stations, stats = genetic_algorithm(mapa, 50, 0.1, num_estacoes)
		elif algoritmo == 'Escalada':
			stations, cost, counter = hill_climbing(mapa, num_estacoes)
			stats = {
				'estados_gerados': counter['estados_gerados'],
				'expansoes': counter['expansoes'],
				'custo_total': 1000 * len(stations) + int(cost * 100),
				'custo_medio': math.floor(cost * 100) / 100.0,
				'tempo_execucao': None  #tempo sera calculado depois 
			}
		elif algoritmo == 'Melhor Primeiro':
			stations, stats = best_first_search(mapa, stations)
		else:
			print("Algoritmo desconhecido.")
			return
		
		if stations is not None:
			stats['tempo_execucao'] = (time.time() - start_time) * 1000
			if stats['custo_medio'] >= 3:
				print(f"\nMatriz {idx + 1} - Nenhuma solução válida encontrada.")
			else:
				print_matrix_and_stats(mapa, stations, stats, idx)
		else:
			if stats is None:
				stats = {
					'estados_gerados': 'N/A',
					'expansoes': 'N/A',
					'custo_total': 'N/A',
					'custo_medio': 'N/A',
					'tempo_execucao': (time.time() - start_time) * 1000
				}
			print(f"Matriz {idx + 1} - Nenhuma solução válida encontrada.")
			print(f"Total de estados gerados: {stats['estados_gerados']}")
			print(f"Total de expansões: {stats['expansoes']}")
			print(f"Custo total: {stats['custo_total']}")
			print(f"Custo médio: {stats['custo_medio']}")
			print(f"Tempo de execução (ms): {stats['tempo_execucao']}")
			
		
			

algoritmo = input("Qual algoritmo quer usar? (A*, Genético, Escalada, Melhor Primeiro): ")
coordenadas_especificas = input("Quer usar coordenadas específicas? (s/n): ").strip().lower() == 's'
num_estacoes_por_matriz = [1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4]  #quantidade de estaçoes por matriz, estao ordenadas.
processar_matrizes(algoritmo, instancias, num_estacoes_por_matriz, coordenadas_especificas)
