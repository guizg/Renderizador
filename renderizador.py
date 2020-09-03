# Desenvolvido por: Luciano Soares <lpsoares@insper.edu.br>
# Disciplina: Computação Gráfica
# Data: 28 de Agosto de 2020

import argparse     # Para tratar os parâmetros da linha de comando
import x3d          # Faz a leitura do arquivo X3D, gera o grafo de cena e faz traversal
import interface    # Janela de visualização baseada no Matplotlib
import gpu          # Simula os recursos de uma GPU

def lineCrossLine(x0, y0, x1, y1, x2, y2, x3, y3):
    if (max(x0,x1) < min(x2,x3)):
        return False 
    
    a1 = (y0-y1)/(x0-x1)
    a2 = (y2-y3)/(x2-x3)
    b1 = y0 - a1*x0
    b2 = y2 - a2*x2

    if a1 == a2:
        return False
    
    xa = (b2 - b1) / (a1 - a2)

    if ( (xa < max( min(x0,x1), min(x2,x3) )) or
        (xa > min( max(x0,x1), max(x2,x3) )) ):
        return False 
    else:
        return True

def shouldPaintPixel(x, y, x0, y0, x1, y1):
    diamondPoints = [(x+0.5, y), (x, y+0.5), (x+1, y+0.5), (x+0.5, y+1)]
    diamondLines = [[diamondPoints[0], diamondPoints[2]], [diamondPoints[2], diamondPoints[3]], [diamondPoints[3], diamondPoints[1]], [diamondPoints[1], diamondPoints[0]]]

    for line in diamondLines:
        if lineCrossLine(x0, y0, x1, y1, line[0][0], line[0][1], line[1][0], line[1][1]):
            return True
    
    return False

def polypoint2D(point, color):
    """ Função usada para renderizar Polypoint2D. """
    color = [color[0]*255, color[1]*255, color[2]*255]

    for i in range(len(point)):
        point[i] = int(point[i])

    for x in range(gpu.GPU.width):
        for y in range(0, gpu.GPU.height):
            x_point_index = 0
            y_point_index = 1

            while x_point_index != len(point):
                if point[x_point_index] == x and point[y_point_index] == y:
                     gpu.GPU.set_pixel(x, y, color[0], color[1], color[2])

                x_point_index += 2
                y_point_index += 2

def polyline2D(lineSegments, color):
    """ Função usada para renderizar Polyline2D. """
    
    color = [color[0]*255, color[1]*255, color[2]*255]

    for x in range(gpu.GPU.width):
        for y in range(gpu.GPU.height):
            if shouldPaintPixel(x, y, lineSegments[0], lineSegments[1], lineSegments[2], lineSegments[3]):
                gpu.GPU.set_pixel(x, y, color[0], color[1], color[2])


# funçoes daqui https://www.geeksforgeeks.org/check-whether-a-given-point-lies-inside-a-triangle-or-not/
def areaTriangle(x1, y1, x2, y2, x3, y3): 
  
    return abs((x1 * ((y2) - y3) + x2 * (y3 - y1)  
                + x3 * (y1 - y2)) / 2.0) 

def isInside(x1, y1, x2, y2, x3, y3, x, y): 
  
    A = areaTriangle(x1, y1, x2, y2, x3, y3) 

    A1 = areaTriangle(x, y, x2, y2, x3, y3) 
      
    A2 = areaTriangle(x1, y1, x, y, x3, y3) 
      
    A3 = areaTriangle(x1, y1, x2, y2, x, y) 
      
    if(round(A,4) == round(A1 + A2 + A3, 4)): 
        return True
    else:
        return False

def triangleSet2D(vertices, color):
    """ Função usada para renderizar TriangleSet2D. """
    color = [color[0]*255, color[1]*255, color[2]*255]

    for x in range(gpu.GPU.width):
        for y in range(gpu.GPU.height):

            if isInside(vertices[0], vertices[1], vertices[2], vertices[3], vertices[4], vertices[5], x+0.5, y+0.5):
                gpu.GPU.set_pixel(x, y, color[0], color[1], color[2])

def triangleSet(point, color):
    """ Função usada para renderizar TriangleSet. """
    # Nessa função você receberá pontos no parâmetro point, esses pontos são uma lista
    # de pontos x, y, e z sempre na ordem. Assim point[0] é o valor da coordenada x do
    # primeiro ponto, point[1] o valor y do primeiro ponto, point[2] o valor z da 
    # coordenada z do primeiro ponto. Já point[3] é a coordenada x do segundo ponto e
    # assim por diante.
    # No TriangleSet os triângulos são informados individualmente, assim os três
    # primeiros pontos definem um triângulo, os três próximos pontos definem um novo
    # triângulo, e assim por diante.
    
    # O print abaixo é só para vocês verificarem o funcionamento, deve ser removido.
    print("TriangleSet : pontos = {0}".format(point)) # imprime no terminal pontos

def viewpoint(position, orientation, fieldOfView):
    """ Função usada para renderizar (na verdade coletar os dados) de Viewpoint. """
    # Na função de viewpoint você receberá a posição, orientação e campo de visão da
    # câmera virtual. Use esses dados para poder calcular e criar a matriz de projeção
    # perspectiva para poder aplicar nos pontos dos objetos geométricos.

    # O print abaixo é só para vocês verificarem o funcionamento, deve ser removido.
    print("Viewpoint : position = {0}, orientation = {1}, fieldOfView = {2}".format(position, orientation, fieldOfView)) # imprime no terminal

def transform(translation, scale, rotation):
    """ Função usada para renderizar (na verdade coletar os dados) de Transform. """
    # A função transform será chamada quando se entrar em um nó X3D do tipo Transform
    # do grafo de cena. Os valores passados são a escala em um vetor [x, y, z]
    # indicando a escala em cada direção, a translação [x, y, z] nas respectivas
    # coordenadas e finalmente a rotação por [x, y, z, t] sendo definida pela rotação
    # do objeto ao redor do eixo x, y, z por t radianos, seguindo a regra da mão direita.
    # Quando se entrar em um nó transform se deverá salvar a matriz de transformação dos
    # modelos do mundo em alguma estrutura de pilha.

    # O print abaixo é só para vocês verificarem o funcionamento, deve ser removido.
    print("Transform : ", end = '')
    if translation:
        print("translation = {0} ".format(translation), end = '') # imprime no terminal
    if scale:
        print("scale = {0} ".format(scale), end = '') # imprime no terminal
    if rotation:
        print("rotation = {0} ".format(rotation), end = '') # imprime no terminal
    print("")

def _transform():
    """ Função usada para renderizar (na verdade coletar os dados) de Transform. """
    # A função _transform será chamada quando se sair em um nó X3D do tipo Transform do
    # grafo de cena. Não são passados valores, porém quando se sai de um nó transform se
    # deverá recuperar a matriz de transformação dos modelos do mundo da estrutura de
    # pilha implementada.

    # O print abaixo é só para vocês verificarem o funcionamento, deve ser removido.
    print("Saindo de Transform")

def triangleStripSet(point, stripCount, color):
    """ Função usada para renderizar TriangleStripSet. """
    # A função triangleStripSet é usada para desenhar tiras de triângulos interconectados,
    # você receberá as coordenadas dos pontos no parâmetro point, esses pontos são uma
    # lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor da coordenada x
    # do primeiro ponto, point[1] o valor y do primeiro ponto, point[2] o valor z da
    # coordenada z do primeiro ponto. Já point[3] é a coordenada x do segundo ponto e assim
    # por diante. No TriangleStripSet a quantidade de vértices a serem usados é informado
    # em uma lista chamada stripCount (perceba que é uma lista).

    # O print abaixo é só para vocês verificarem o funcionamento, deve ser removido.
    print("TriangleStripSet : pontos = {0} ".format(point), end = '') # imprime no terminal pontos
    for i, strip in enumerate(stripCount):
        print("strip[{0}] = {1} ".format(i, strip), end = '') # imprime no terminal
    print("")

def indexedTriangleStripSet(point, index, color):
    """ Função usada para renderizar IndexedTriangleStripSet. """
    # A função indexedTriangleStripSet é usada para desenhar tiras de triângulos
    # interconectados, você receberá as coordenadas dos pontos no parâmetro point, esses
    # pontos são uma lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor
    # da coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto, point[2]
    # o valor z da coordenada z do primeiro ponto. Já point[3] é a coordenada x do
    # segundo ponto e assim por diante. No IndexedTriangleStripSet uma lista informando
    # como conectar os vértices é informada em index, o valor -1 indica que a lista
    # acabou. A ordem de conexão será de 3 em 3 pulando um índice. Por exemplo: o
    # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
    # depois 2, 3 e 4, e assim por diante.
    
    # O print abaixo é só para vocês verificarem o funcionamento, deve ser removido.
    print("IndexedTriangleStripSet : pontos = {0}, index = {1}".format(point, index)) # imprime no terminal pontos

def box(size, color):
    """ Função usada para renderizar Boxes. """
    # A função box é usada para desenhar paralelepípedos na cena. O Box é centrada no
    # (0, 0, 0) no sistema de coordenadas local e alinhado com os eixos de coordenadas
    # locais. O argumento size especifica as extensões da caixa ao longo dos eixos X, Y
    # e Z, respectivamente, e cada valor do tamanho deve ser maior que zero. Para desenha
    # essa caixa você vai provavelmente querer tesselar ela em triângulos, para isso
    # encontre os vértices e defina os triângulos.

    # O print abaixo é só para vocês verificarem o funcionamento, deve ser removido.
    print("Box : size = {0}".format(size)) # imprime no terminal pontos

def indexedFaceSet(coord, coordIndex, colorPerVertex, color, colorIndex, texCoord, texCoordIndex, current_color, current_texture):
    """ Função usada para renderizar IndexedFaceSet. """
    # A função indexedFaceSet é usada para desenhar malhas de triângulos. Ela funciona de
    # forma muito simular a IndexedTriangleStripSet porém com mais recursos.
    # Você receberá as coordenadas dos pontos no parâmetro cord, esses
    # pontos são uma lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor
    # da coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto, point[2]
    # o valor z da coordenada z do primeiro ponto. Já point[3] é a coordenada x do
    # segundo ponto e assim por diante. No IndexedFaceSet uma lista informando
    # como conectar os vértices é informada em coordIndex, o valor -1 indica que a lista
    # acabou. A ordem de conexão será de 3 em 3 pulando um índice. Por exemplo: o
    # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
    # depois 2, 3 e 4, e assim por diante.
    # Adicionalmente essa implementação do IndexedFace suport cores por vértices, assim
    # a se a flag colorPerVertex estiver habilidades, os vértices também possuirão cores
    # que servem para definir a cor interna dos poligonos, para isso faça um cálculo
    # baricêntrico de que cor deverá ter aquela posição. Da mesma forma se pode definir uma
    # textura para o poligono, para isso, use as coordenadas de textura e depois aplique a
    # cor da textura conforme a posição do mapeamento. Dentro da classe GPU já está
    # implementadado um método para a leitura de imagens.
    
    # O print abaixo é só para vocês verificarem o funcionamento, deve ser removido.
    print("IndexedFaceSet : ")
    if coord:
        print("\tpontos(x, y, z) = {0}, coordIndex = {1}".format(coord, coordIndex)) # imprime no terminal
    if colorPerVertex:
        print("\tcores(r, g, b) = {0}, colorIndex = {1}".format(color, colorIndex)) # imprime no terminal
    if texCoord:
        print("\tpontos(u, v) = {0}, texCoordIndex = {1}".format(texCoord, texCoordIndex)) # imprime no terminal
    if(current_texture):
        image = gpu.GPU.load_texture(current_texture[0])
        print("\t Matriz com image = {0}".format(image))

# Defina o tamanhã da tela que melhor sirva para perceber a renderização
LARGURA = 30
ALTURA = 20

if __name__ == '__main__':

    # Valores padrão da aplicação
    width = LARGURA
    height = ALTURA
<<<<<<< HEAD
    x3d_file = "exemplo9.x3d"
=======
<<<<<<< HEAD
    x3d_file = "exemplo4.x3d"
=======
    x3d_file = "exemplo3.x3d"
>>>>>>> p1 done
>>>>>>> p1 done
    image_file = "tela.png"

    # Tratando entrada de parâmetro
    parser = argparse.ArgumentParser(add_help=False)   # parser para linha de comando
    parser.add_argument("-i", "--input", help="arquivo X3D de entrada")
    parser.add_argument("-o", "--output", help="arquivo 2D de saída (imagem)")
    parser.add_argument("-w", "--width", help="resolução horizonta", type=int)
    parser.add_argument("-h", "--height", help="resolução vertical", type=int)
    parser.add_argument("-q", "--quiet", help="não exibe janela de visualização", action='store_true')
    args = parser.parse_args() # parse the arguments
    if args.input: x3d_file = args.input
    if args.output: image_file = args.output
    if args.width: width = args.width
    if args.height: height = args.height

    # Iniciando simulação de GPU
    gpu.GPU(width, height, image_file)

    # Abre arquivo X3D
    scene = x3d.X3D(x3d_file)
    scene.set_resolution(width, height)

    # funções que irão fazer o rendering
    x3d.X3D.render["Polypoint2D"] = polypoint2D
    x3d.X3D.render["Polyline2D"] = polyline2D
    x3d.X3D.render["TriangleSet2D"] = triangleSet2D
    x3d.X3D.render["TriangleSet"] = triangleSet
    x3d.X3D.render["Viewpoint"] = viewpoint
    x3d.X3D.render["Transform"] = transform
    x3d.X3D.render["_Transform"] = _transform
    x3d.X3D.render["TriangleStripSet"] = triangleStripSet
    x3d.X3D.render["IndexedTriangleStripSet"] = indexedTriangleStripSet
    x3d.X3D.render["Box"] = box
    x3d.X3D.render["IndexedFaceSet"] = indexedFaceSet

    # Se no modo silencioso não configurar janela de visualização
    if not args.quiet:
        window = interface.Interface(width, height)
        scene.set_preview(window)

    scene.parse() # faz o traversal no grafo de cena

    # Se no modo silencioso salvar imagem e não mostrar janela de visualização
    if args.quiet:
        gpu.GPU.save_image() # Salva imagem em arquivo
    else:
        window.image_saver = gpu.GPU.save_image # pasa a função para salvar imagens
        window.preview(gpu.GPU._frame_buffer) # mostra janela de visualização
