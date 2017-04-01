%% PREGUNTA 1
%%
%Leer fotos del fondo

sdirectory = 'fondo';
files = dir([sdirectory '/*.jpg']);
N  = length(files);
fondo = cell(1,N);
for k=1:N
    filename = [sdirectory '/' files(k).name];
    Irgb = imread(filename);
    Igray = rgb2gray(Irgb);
    fondo{k} = Igray;
    
end

%% Leer seq1


s1directory = 'seq1';
files1 = dir([s1directory '/*.jpg']);
N  = length(files1);
s1q = cell(1,N);
for k=1:N
    filename = [s1directory '/' files1(k).name];
    Irgb = imread(filename);
    Igray = rgb2gray(Irgb);
    s1q{k} = Igray;
    
end

%% Leer seq2
s1directory = 'seq2';
files2 = dir([s1directory '/*.jpg']);
N  = length(files2);
s2q = cell(1,N);
for k=1:N
    filename = [s1directory '/' files2(k).name];
    Irgb = imread(filename);
    Igray = rgb2gray(Irgb);
    s2q{k} = Igray;
    
end

%% Leer seq3

s1directory = 'seq3';
files3 = dir([s1directory '/*.jpg']);
N  = length(files3);
s3q = cell(1,N);
for k=1:N
    filename = [s1directory '/' files3(k).name];
    Irgb = imread(filename);
    Igray = rgb2gray(Irgb);
    s3q{k} = Igray;
    
end

%% Queremos hacer la diferencia ahora

%Primera secuencia
v = VideoWriter('s1q.avi');
v.FrameRate = 20;
tres= 30;
S = cell(1,499);
open(v);
for k=2:500
    delta=s1q{k} - s1q{k-1} > tres;
    S{k-1} = delta;
    writeVideo(v,uint8(S{k-1})*255);
 
end 
close(v);

%%

% Si le achicas el umbral implica mas ruido
%Segunda secuencia
w = VideoWriter('s2q.avi');
w.FrameRate = 20;
tres2= 30;
S2 = cell(1,499);
open(w);
for k=2:400
    delta1=s2q{k} - s2q{k-1} > tres2;
    S2{k-1} = delta1;
    writeVideo(w,uint8(S2{k-1})*255);
 
end 
close(w);

%%
%Tercera secuencia
z = VideoWriter('s3q.avi');
z.FrameRate = 20;
tres3= 30;
S3 = cell(1,174);
open(z);
for k=2:175
    delta2=s3q{k} - s3q{k-1} > tres3;
    S3{k-1} = delta2;
    writeVideo(z,uint8(S3{k-1}*255));
 
end 
close(z);
 
%% PREGUNTA 2
%% Promedio
prom = zeros(800,1280);
for k=1:length(files)
    prom = prom + double((fondo{k}));
end
promedio=prom./175;
prom2 = uint8(prom./175);
%imshow(prom2)
%% Desviacion estandar
desv = zeros(800,1280);
for k = 1:length(files)
    desv = desv + ((double(fondo{k}) - promedio).^2) ;
end
desv = desv./174;
desv = sqrt(desv);
%imshow(uint8(desv));
%%
alpha = 20;
limSuperior = promedio + double(alpha)*desv;
limInferior = promedio- double(alpha)*desv;
%%
deteccion1 = cell(1,499);
d1 = VideoWriter('deteccion1.avi');
d1.FrameRate = 20;

open(d1);
for i=1:length(files1)
    deteccion1{1,i} = (s1q{i} > limSuperior | s1q{i} < limInferior );
    writeVideo(d1,uint8(deteccion1{i})*255);
end
close(d1);
%%
deteccion2 = cell(1,400);
d2 = VideoWriter('deteccion2.avi');
d2.FrameRate = 20;

open(d2);
for i=1:length(files2)
    deteccion2{i} = (s2q{i} > limSuperior | s2q{i} < limInferior );
    writeVideo(d2,uint8(deteccion2{i})*255);
end
close(d2);

%%
deteccion3 = cell(1,174);
d3 = VideoWriter('deteccion3.avi');
d3.FrameRate = 20;

open(d3);
for i=1:length(files3)
    deteccion3{i} = (s3q{i} > limSuperior | s3q{i} < limInferior );
    writeVideo(d3,uint8(deteccion3{i})*255);
end
close(d3);


%%
%PREGUNTA 3 
%% Sequencia 1 Filas
v=cell(1,500);
for j=1:500
    v{1,j} = zeros(800,1);
end

for z=1:length(deteccion1)
    for i=1:800
        for j=1:1280
            v{1,z}(i,1) = v{1,z}(i,1) + deteccion1{z}(i,j);
        end
    end
end

%% Seq 1 Columnas
vcol = cell(1,500);
for j=1:500
    vcol{1,j} = zeros(1280,1);
end

for z=1:length(deteccion1)
    for j=1:1280
        for i=1:800
            vcol{1,z}(j,1) = vcol{1,z}(j,1) + deteccion1{z}(i,j);
        end
    end
end

%% Histograma seq1 Filas
cant1 = 1:800;
for i=1:500
    plot(cant1,v{i})
    ylim([0,150]);
    %pause
    %disp(i)
end
%% Histograma seq1 Columnas
cant11 = 1:1280;
for i=1:500
    plot(cant11,vcol{i})
    ylim([0,150]);
%     pause
%     disp(i)
end

%% Calculo de boundaries de fila para cada frame seq 1
y1Rows = zeros(1,500);
y2Rows = zeros(1,500);
for fr=1:500
    [a,b] = max(v{fr});
    linea = ones(1,100)*b;
    y = linspace(0,a);
    plot(cant1,v{fr})
    hold on
    plot(linea,y);
    beta =50;
    ileft=b;
    iright = b;
    for i=1:length(cant1)
        if(i<b) 
            if abs(v{fr}(b,1) - v{fr}(i,1)) > beta
                ileft = i;
            end
        else
            if abs(v{fr}(b,1) - v{fr}(i,1))> beta
                iright = i;
                break
            end
        end
    end 
    y1Rows(fr) = ileft;
    y2Rows(fr) = iright;
    lineaLeft = ones(1,100)*ileft;
    yLeft = linspace(0,v{fr}(ileft,1));
    plot(lineaLeft,yLeft,'r');
    lineaRight = ones(1,100)*iright;
    yRight = linspace(0,v{fr}(iright,1));
    plot(lineaRight,yRight,'r')
    %pause
    hold off
    
end
%% %% Calculo de boundaries de columna para cada frame seq 1
x1Col = zeros(1,500);
x2Col = zeros(1,500);
for fr=1:500
    [cCols1,dCols1] = max(vcol{fr});
    lineaCols1= ones(1,100)*dCols1;
    yCols1 = linspace(0,cCols1);
    plot(cant11,vcol{fr})
    hold on
    plot(lineaCols1,yCols1);
    betaCols1 = 50;
    ileftCols1=dCols1;
    irightCols1 = dCols1;
    for i=1:length(cant11)
        if(i<dCols1) 
            if abs(vcol{fr}(dCols1,1) - vcol{fr}(i,1)) > betaCols1
                ileftCols1 = i;
            end
        else
            if abs(vcol{fr}(dCols1,1) - vcol{fr}(i,1))> betaCols1
                irightCols1 = i;
                break
            end
        end
    end 
    x1Col(fr) = ileftCols1;
    x2Col(fr) = irightCols1;
    lineaLeftCols1 = ones(1,100)*ileftCols1;
    yLeftCols1 = linspace(0,vcol{fr}(ileftCols1,1));
    plot(lineaLeftCols1,yLeftCols1,'r');
    lineaRightCols1 = ones(1,100)*irightCols1;
    yRightCols1 = linspace(0,vcol{fr}(irightCols1,1));
    plot(lineaRightCols1,yRightCols1,'r')
    %pause
    hold off
    
end

%% SEQ1 BLOB- Creacion de blob de deteccion para primera secuencia 

caja1avi = VideoWriter('caja1.avi');
caja1avi.FrameRate = 20;
open(caja1avi);
for i=1:500

    vertices = [x1Col(i),y1Rows(i), y2Rows(i)- y1Rows(i), x2Col(i)- x1Col(i)  ];
    caja1 = uint8(255*deteccion1{1,i});
    imgDet = insertShape(caja1,'rectangle',vertices);
    writeVideo(caja1avi,uint8(imgDet)*255);
end
close(caja1avi);
%%
%% Parte 4 Seq 2
Estevez= 0.5;
%Tengo xiCol3, yiRows3
Ancho1 = x2Col - x1Col;
Largo1 = y2Rows - y1Rows;
%Definimos los nuevos vectores
PredictX1 = zeros(1,500); 
PredictY1 = zeros(1,500);
PredictW1 = zeros(1,500);
PredictH1 = zeros(1,500);
PredictX1(1) = x1Col(1);
PredictX1(2) = x1Col(2);
PredictY1(1) = y1Rows(1);
PredictY1(2) = y1Rows(2);
PredictW1(1) = Ancho1(1);
PredictW1(2) = Ancho1(2);
PredictH1(1) = Largo1(2);
PredictH1(2) = Largo1(2);
for i=3:500
    PredictX1(i) = x1Col(i-1) + Estevez*(   x1Col(i-1) - x1Col(i-2));
    PredictY1(i) = y1Rows(i-1) + Estevez*(   y1Rows(i-1) - y1Rows(i-2));
    PredictW1(i) = Ancho1(i-1) + Estevez*(   Ancho1(i-1) - Ancho1(i-2));
    PredictH1(i) = Largo1(i-1) + Estevez*(   Largo1(i-1) - Largo1(i-2));
end
%%
%% Parte 4 SEQ2 BLOB- TRACKING

caja1Trackavi = VideoWriter('caja1Track.avi');
caja1Trackavi.FrameRate = 20;
open(caja1Trackavi);
for i=1:500
    vertices1 = [x1Col(i),y1Rows(i), x2Col(i)- x1Col(i),y2Rows(i)- y1Rows(i)  ];
    verticesTrack1 = [PredictX1(i),PredictY1(i), PredictW1(i),PredictH1(i)  ];
    caja1 = uint8(255*deteccion1{i});
%     imgDetTrack = insertShape(caja3Track,'rectangle',{vertices, verticesTrack}, 'Color','blue');
    img1 = insertShape(caja1,'rectangle', vertices1);
    imgDetTrack1 = insertShape(img1,'rectangle', verticesTrack1, 'Color', 'blue');
    writeVideo(caja1Trackavi,uint8(imgDetTrack1)*255);
end
close(caja1Trackavi);


%%
%%SEQUENCE 2
%% Sequencia 2 Filas
vSeq2=cell(1,400);
for j=1:400
    vSeq2{1,j} = zeros(800,1);
end

for z=1:length(deteccion2)
    for i=1:800
        for j=1:1280
            vSeq2{1,z}(i,1) = vSeq2{1,z}(i,1) + deteccion2{z}(i,j);
        end
    end
end
%% Mucho Ruido Seq2, cambiar!!!
%% Seq 2 Columnas
vcolSeq2 = cell(1,400);
for j=1:400
    vcolSeq2{1,j} = zeros(1280,1);
end

for z=1:length(deteccion2)
    for j=1:1280
        for i=1:800
            vcolSeq2{1,z}(j,1) = vcolSeq2{1,z}(j,1) + deteccion2{z}(i,j);
        end
    end
end
%% Histograma seq2 Filas
cant2 = 1:800;
for i=1:400
    plot(cant2,vSeq2{i})
    ylim([0,250]);
%     pause
%     disp(i)
end
%% Histograma seq2 Columnas
cant21 = 1:1280;
for i=1:400
    plot(cant21,vcolSeq2{i})
    ylim([0,150]);
    %pause
    %disp(i)
end
%% Calculo de boundaries de fila para cada frame seq 2
y1Rows2 = zeros(1,400);
y2Rows2 = zeros(1,400);
for fr=1:400
    [aRows2,bRows2] = max(vSeq2{fr});
    lineaRows2 = ones(1,100)*bRows2;
    yRows2 = linspace(0,aRows2);
    plot(cant2,vSeq2{fr})
    hold on
    plot(lineaRows2,yRows2);
    betaRows2 =100;
    ileftRows2=bRows2;
    irightRows2 = bRows2;
    for i=1:400
        if(i<bRows2) 
            if abs(vSeq2{fr}(bRows2,1) - vSeq2{fr}(i,1)) > betaRows2
                ileftRows2 = i;
            end
        else
            if abs(vSeq2{fr}(bRows2,1) - vSeq2{fr}(i,1))> betaRows2
                irightRows2 = i;
                break
            end
        end
    end 
    y1Rows2(fr) = ileftRows2;
    y2Rows2(fr) = irightRows2;
    lineaLeftRows2 = ones(1,100)*ileftRows2;
    yLeftRows2 = linspace(0,vSeq2{fr}(ileftRows2,1));
    plot(lineaLeftRows2,yLeftRows2,'r');
    lineaRightRows2 = ones(1,100)*irightRows2;
    yRightRows2 = linspace(0,vSeq2{fr}(irightRows2,1));
    plot(lineaRightRows2,yRightRows2,'r')
    %pause
    hold off
    
end
%% %% Calculo de boundaries de columna para cada frame seq 2
x1Col2 = zeros(1,400);
x2Col2 = zeros(1,400);
for fr=1:400
    [cCols2,dCols2] = max(vcolSeq2{fr});
    lineaCols2= ones(1,100)*dCols2;
    yCols2 = linspace(0,cCols2);
    plot(cant21,vcolSeq2{fr})
    hold on
    plot(lineaCols2,yCols2);
    betaCols2 = 100;
    ileftCols2=dCols2;
    irightCols2 = dCols2;
    for i=1:length(cant21)
        if(i<dCols2) 
            if abs(vcolSeq2{fr}(dCols2,1) - vcolSeq2{fr}(i,1)) > betaCols2
                ileftCols2 = i;
            end
        else
            if abs(vcolSeq2{fr}(dCols2,1) - vcolSeq2{fr}(i,1))> betaCols2
                irightCols2 = i;
                break
            end
        end
    end 
    x1Col2(fr) = ileftCols2;
    x2Col2(fr) = irightCols2;
    lineaLeftCols2 = ones(1,100)*ileftCols2;
    yLeftCols2 = linspace(0,vcolSeq2{fr}(ileftCols2,1));
    plot(lineaLeftCols2,yLeftCols2,'r');
    lineaRightCols2 = ones(1,100)*irightCols2;
    yRightCols2 = linspace(0,vcolSeq2{fr}(irightCols2,1));
    plot(lineaRightCols2,yRightCols2,'r')
    %pause
    hold off
    
end

%% SEQ2 BLOB- Creacion de blob de deteccion para segunda secuencia 

caja2avi = VideoWriter('caja2.avi');
caja2avi.FrameRate = 20;
open(caja2avi);
for i=1:400

    vertices = [x1Col2(i),y1Rows2(i), x2Col2(i)- x1Col2(i),y2Rows2(i)- y1Rows2(i)  ];
    caja2 = uint8(255*deteccion2{1,i});
    imgDet = insertShape(caja2,'rectangle',vertices);
    writeVideo(caja2avi,uint8(imgDet)*255);
end
close(caja2avi);
%% Parte 4 Seq 2
Estevez= 0.5;
%Tengo xiCol3, yiRows3
Ancho2 = x2Col2 - x1Col2;
Largo2 = y2Rows2 - y1Rows2;
%Definimos los nuevos vectores
PredictX2 = zeros(1,400); 
PredictY2 = zeros(1,400);
PredictW2 = zeros(1,400);
PredictH2 = zeros(1,400);
PredictX2(1) = x1Col2(1);
PredictX2(2) = x1Col2(2);
PredictY2(1) = y1Rows2(1);
PredictY2(2) = y1Rows2(2);
PredictW2(1) = Ancho2(1);
PredictW2(2) = Ancho2(2);
PredictH2(1) = Largo2(2);
PredictH2(2) = Largo2(2);
for i=3:400
    PredictX2(i) = x1Col2(i-1) + Estevez*(   x1Col2(i-1) - x1Col2(i-2));
    PredictY2(i) = y1Rows2(i-1) + Estevez*(   y1Rows2(i-1) - y1Rows2(i-2));
    PredictW2(i) = Ancho2(i-1) + Estevez*(   Ancho2(i-1) - Ancho2(i-2));
    PredictH2(i) = Largo2(i-1) + Estevez*(   Largo2(i-1) - Largo2(i-2));
end
%%
%% Parte 4 SEQ2 BLOB- TRACKING

caja2Trackavi = VideoWriter('caja2Track.avi');
caja2Trackavi.FrameRate = 20;
open(caja2Trackavi);
for i=1:400
    vertices2 = [x1Col2(i),y1Rows2(i), x2Col2(i)- x1Col2(i),y2Rows2(i)- y1Rows2(i)  ];
    verticesTrack2 = [PredictX2(i),PredictY2(i), PredictW2(i),PredictH2(i)  ];
    caja2 = uint8(255*deteccion2{i});
%     imgDetTrack = insertShape(caja3Track,'rectangle',{vertices, verticesTrack}, 'Color','blue');
    img2 = insertShape(caja2,'rectangle', vertices2);
    imgDetTrack2 = insertShape(img2,'rectangle', verticesTrack2, 'Color', 'blue');
    writeVideo(caja2Trackavi,uint8(imgDetTrack2)*255);
end
close(caja2Trackavi);

%%
%% SEQUENCE 3
%% Sequencia 3 Filas
vSeq3=cell(1,175);
for j=1:175
    vSeq3{1,j} = zeros(800,1);
end

for z=1:length(deteccion3)
    for i=1:800
        for j=1:1280
            vSeq3{1,z}(i,1) = vSeq3{1,z}(i,1) + deteccion3{z}(i,j);
        end
    end
end

%% Seq 3 Columnas
vcolSeq3 = cell(1,175);
for j=1:175
    vcolSeq3{1,j} = zeros(1280,1);
end

for z=1:length(deteccion3)
    for j=1:1280
        for i=1:800
            vcolSeq3{1,z}(j,1) = vcolSeq3{1,z}(j,1) + deteccion3{z}(i,j);
        end
    end
end
%% Histograma seq3 Filas
cant3 = 1:800;
for i=1:175
    plot(cant3,vSeq3{i})
    ylim([0,250]);
    %pause
    %disp(i)
end
%% Histograma seq3 Columnas
cant31 = 1:1280;
for i=1:175
    plot(cant31,vcolSeq3{i})
    ylim([0,150]);
%     pause
%     disp(i)
end
%% Calculo de boundaries de fila para cada frame seq 3
y1Rows3 = zeros(1,175);
y2Rows3 = zeros(1,175);
for fr=1:175
    [aRows3,bRows3] = max(vSeq3{fr});
    lineaRows3 = ones(1,100)*bRows3;
    yRows3 = linspace(0,aRows3);
    plot(cant3,vSeq3{fr})
    hold on
    plot(lineaRows3,yRows3);
    betaRows3 =35;
    ileftRows3=bRows3;
    irightRows3 = bRows3;
    for i=1:length(cant3)
        if(i<bRows3) 
            if abs(vSeq3{fr}(bRows3,1) - vSeq3{fr}(i,1)) > betaRows3
                ileftRows3 = i;
            end
        else
            if abs(vSeq3{fr}(bRows3,1) - vSeq3{fr}(i,1))> betaRows3
                irightRows3 = i;
                break
            end
        end
    end 
    y1Rows3(fr) = ileftRows3;
    y2Rows3(fr) = irightRows3;
    lineaLeftRows3 = ones(1,100)*ileftRows3;
    yLeftRows3 = linspace(0,vSeq3{fr}(ileftRows3,1));
    plot(lineaLeftRows3,yLeftRows3,'r');
    lineaRightRows3 = ones(1,100)*irightRows3;
    yRightRows3 = linspace(0,vSeq3{fr}(irightRows3,1));
    plot(lineaRightRows3,yRightRows3,'r')
    %pause
    hold off
    
end
%% %% Calculo de boundaries de columna para cada frame seq 3
x1Col3 = zeros(1,175);
x2Col3 = zeros(1,175);
for fr=1:175
    [cCols3,dCols3] = max(vcolSeq3{fr});
    lineaCols3= ones(1,100)*dCols3;
    yCols3 = linspace(0,cCols3);
    plot(cant31,vcolSeq3{fr})
    hold on
    plot(lineaCols3,yCols3);
    betaCols3 = 25;
    ileftCols3=dCols3;
    irightCols3 = dCols3;
    for i=1:length(cant31)
        if(i<dCols3) 
            if abs(vcolSeq3{fr}(dCols3,1) - vcolSeq3{fr}(i,1)) > betaCols3
                ileftCols3 = i;
            end
        else
            if abs(vcolSeq3{fr}(dCols3,1) - vcolSeq3{fr}(i,1))> betaCols3
                irightCols3 = i;
                break
            end
        end
    end 
    x1Col3(fr) = ileftCols3;
    x2Col3(fr) = irightCols3;
    lineaLeftCols3 = ones(1,100)*ileftCols3;
    yLeftCols3 = linspace(0,vcolSeq3{fr}(ileftCols3,1));
    plot(lineaLeftCols3,yLeftCols3,'r');
    lineaRightCols3 = ones(1,100)*irightCols3;
    yRightCols3 = linspace(0,vcolSeq3{fr}(irightCols3,1));
    plot(lineaRightCols3,yRightCols3,'r')
    %pause
    hold off
    
end
%%
%% SEQ3 BLOB- Creacion de blob de deteccion para tercera secuencia 

caja3avi = VideoWriter('caja3.avi');
caja3avi.FrameRate = 12;
open(caja3avi);
for i=1:175

    vertices = [x1Col3(i),y1Rows3(i), x2Col3(i)- x1Col3(i),y2Rows3(i)- y1Rows3(i)  ];
    caja3 = uint8(255*deteccion3{1,i});
    imgDet = insertShape(caja3,'rectangle',vertices);
    writeVideo(caja3avi,uint8(imgDet)*255);
end
close(caja3avi);
%% Parte 4 Seq 3
Estevez= 0.5;
%Tengo xiCol3, yiRows3
Ancho3 = x2Col3 - x1Col3;
Largo3 = y2Rows3 - y1Rows3;
%Definimos los nuevos vectores
PredictX3 = zeros(1,175); 
PredictY3 = zeros(1,175);
PredictW3 = zeros(1,175);
PredictH3 = zeros(1,175);
PredictX3(1) = x1Col3(1);
PredictX3(2) = x1Col3(2);
PredictY3(1) = y1Rows3(1);
PredictY3(2) = y1Rows3(2);
PredictW3(1) = Ancho3(1);
PredictW3(2) = Ancho3(2);
PredictH3(1) = Largo3(2);
PredictH3(2) = Largo3(2);
for i=3:175
    PredictX3(i) = x1Col3(i-1) + Estevez*(   x1Col3(i-1) - x1Col3(i-2));
    PredictY3(i) = y1Rows3(i-1) + Estevez*(   y1Rows3(i-1) - y1Rows3(i-2));
    PredictW3(i) = Ancho3(i-1) + Estevez*(   Ancho3(i-1) - Ancho3(i-2));
    PredictH3(i) = Largo3(i-1) + Estevez*(   Largo3(i-1) - Largo3(i-2));
end

%% Parte 4 SEQ3 BLOB- TRACKING

caja3Trackavi = VideoWriter('caja3Track.avi');
caja3Trackavi.FrameRate = 12;
open(caja3Trackavi);
for i=1:175
    vertices = [x1Col3(i),y1Rows3(i), x2Col3(i)- x1Col3(i),y2Rows3(i)- y1Rows3(i)  ];
    verticesTrack = [PredictX3(i),PredictY3(i), PredictW3(i),PredictH3(i)  ];
    caja3 = uint8(255*deteccion3{i});
%     imgDetTrack = insertShape(caja3Track,'rectangle',{vertices, verticesTrack}, 'Color','blue');
    img = insertShape(caja3,'rectangle', vertices);
    imgDetTrack = insertShape(img,'rectangle', verticesTrack, 'Color', 'blue');
    writeVideo(caja3Trackavi,uint8(imgDetTrack)*255);
end
close(caja3Trackavi);
%%
