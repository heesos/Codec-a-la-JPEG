close all;
clear;
%% ZMIENNE
I=imread('Tiger.jpg');
I=rgb2ycbcr(I);
[width,height,x]=size(I);

BLOCKSIZE=8;
blockR=zeros([BLOCKSIZE BLOCKSIZE floor(width*height/BLOCKSIZE^2)]);
blockG=zeros([BLOCKSIZE BLOCKSIZE floor(width*height/BLOCKSIZE^2)]);
blockB=zeros([BLOCKSIZE BLOCKSIZE floor(width*height/BLOCKSIZE^2)]);
y=zeros([BLOCKSIZE BLOCKSIZE floor(width*height/BLOCKSIZE^2)]);
cb=zeros([BLOCKSIZE BLOCKSIZE floor(width*height/BLOCKSIZE^2)]);
cr=zeros([BLOCKSIZE BLOCKSIZE floor(width*height/BLOCKSIZE^2)]);

tab_kwantyzacji = [16 11 10 16 24 40 51 61;
                   12 12 14 19 26 58 60 55;
                   14 13 16 24 40 57 69 56;
                   14 17 22 29 51 87 80 62;
                   18 22 37 56 68 109 103 77;
                   24 35 55 64 81 104 113 92;
                   49 64 78 87 103 121 120 101;
                   72 92 95 98 112 100 103 99];
%% PODZIAL NA BLOKI
    n=1;
    for j=1:BLOCKSIZE:height
        if j+BLOCKSIZE>height
             for i = 1:BLOCKSIZE:width
                   if i+BLOCKSIZE>width
                      blockR(:,:,n)=I(i:width,j:height,1);
                      blockG(:,:,n)=I(i:width,j:height,2);
                      blockB(:,:,n)=I(i:width,j:height,3);           
                       break;
                   else 
                      blockR(:,:,n)=I(i:BLOCKSIZE+i-1,j:height,1);
                      blockG(:,:,n)=I(i:BLOCKSIZE+i-1,j:height,2);
                      blockB(:,:,n)=I(i:BLOCKSIZE+i-1,j:height,3);
                   end
                   n=n+1;
             end
             break;
        else
            for i = 1:BLOCKSIZE:width
                   if i+BLOCKSIZE>width
                      blockR(:,:,n)=I(i:width,j:BLOCKSIZE+j-1,1);
                      blockG(:,:,n)=I(i:width,j:BLOCKSIZE+j-1,2);
                      blockB(:,:,n)=I(i:width,j:BLOCKSIZE+j-1,3);
                       n=n+1;
                       break;
                   else 
                       blockR(:,:,n)=I(i:BLOCKSIZE+i-1,j:BLOCKSIZE+j-1,1);
                       blockG(:,:,n)=I(i:BLOCKSIZE+i-1,j:BLOCKSIZE+j-1,2);
                       blockB(:,:,n)=I(i:BLOCKSIZE+i-1,j:BLOCKSIZE+j-1,3);
                       n=n+1;
                   end
           end
        end
    end

numberofBlocks=size(blockR);
%% DCT
for l=1:1:numberofBlocks(3)
    y(:,:,l)=dct2(blockR(:,:,l));
    cb(:,:,l)=dct2(blockG(:,:,l));
    cr(:,:,l)=dct2(blockB(:,:,l));
end
%% KWANTYZACJA 
for i = 1:1:numberofBlocks(3)
    for j = 1:1:BLOCKSIZE
        for k = 1:1:BLOCKSIZE
            y(k,j,i) =  round(y(k,j,i)/tab_kwantyzacji(k,j));
            cr(k,j,i) = round(cr(k,j,i)/tab_kwantyzacji(k,j));
            cb(k,j,i) = round(cb(k,j,i)/tab_kwantyzacji(k,j));             
        end
    end
end


%% ZAMIANIA NA WEKTORY ZIG-ZAG
y3 = zeros([BLOCKSIZE^2 width*height/BLOCKSIZE^2]);
cb3 =zeros([BLOCKSIZE^2 width*height/BLOCKSIZE^2]);
cr3 = zeros([BLOCKSIZE^2 width*height/BLOCKSIZE^2]);
% y3=[];
% cb3=[];
% cr3=[];
for j = 1:1:numberofBlocks(3)
    y3(:,j) = zigzag(y(:,:,j));
    cb3(:,j) = zigzag(cb(:,:,j));
    cr3(:,j) = zigzag(cr(:,:,j));
end

%% HUFFMAN CODING

B =[];
R =[];
G =[];
for i =1:1:numberofBlocks(3)
    for k=1:1:BLOCKSIZE^2
      R(end+1)=y3(k,i);
      G(end+1)=cb3(k,i);
      B(end+1)=cr3(k,i);
   end
end
long = [R,G,B];

symbols = unique(long);
counts = hist(long, symbols);
p = double(counts) ./ sum(counts);
[dict,avglen] = huffmandict(symbols,p); 
comp = huffmanenco(long,dict);


%% ZAPIS DO PLIKU

save('save.jmm','comp','dict');

%% ODCZYT Z PLIKU

load('save.jmm','-mat');
    
%% HUFFMAN DECODING
dhsig1 = huffmandeco(comp,dict);
y5 = dhsig1(1:BLOCKSIZE^2*numberofBlocks(3));
cb5 = dhsig1(1+BLOCKSIZE^2*numberofBlocks(3):2*BLOCKSIZE^2*numberofBlocks(3));
cr5 = dhsig1(2*BLOCKSIZE^2*numberofBlocks(3)+1:3*BLOCKSIZE^2*numberofBlocks(3));


y3=reshape(y5,BLOCKSIZE^2,numberofBlocks(3));
cb3=reshape(cb5,BLOCKSIZE^2,numberofBlocks(3));
cr3=reshape(cr5,BLOCKSIZE^2,numberofBlocks(3));
       

%% ZIG-ZAG ALE OD TYLU 

for i = 1:1:numberofBlocks(3)
    x=y3(:,i).';
    y(:,:,i)=izigsc(x,BLOCKSIZE);
    x1=cb3(:,i).';
    cb(:,:,i)=izigsc(x1,BLOCKSIZE);
    x2=cr3(:,i).';
    cr(:,:,i)=izigsc(x2,BLOCKSIZE);
end

%% Dekompresja
for i = 1:1:numberofBlocks(3)
    for j = 1:1:BLOCKSIZE
        for k = 1:1:BLOCKSIZE
            y(k,j,i) = round(y(k,j,i)*tab_kwantyzacji(k,j));
            cr(k,j,i) = round(cr(k,j,i)*tab_kwantyzacji(k,j));
            cb(k,j,i) = round(cb(k,j,i)*tab_kwantyzacji(k,j));        
        end
    end
end

for p = 1:1:numberofBlocks(3)
    y(:,:,p) = idct2(y(:,:,p));
    cb(:,:,p) = idct2(cb(:,:,p));
    cr(:,:,p) = idct2(cr(:,:,p));
end


%% Łączenie bloków w obraz
decodedR=[];
decodedG=[];
decodedB=[];
numberofColumns=floor(width/BLOCKSIZE);
numberofRows=floor(height/BLOCKSIZE);

    counter=1;
    for i=0:1:numberofRows-1
        for j=0:1:numberofColumns-1      
            decodedR(j*BLOCKSIZE+1:(j+1)*BLOCKSIZE,1+i*BLOCKSIZE:(i+1)*BLOCKSIZE)=y(:,:,counter);
            decodedG(j*BLOCKSIZE+1:(j+1)*BLOCKSIZE,1+i*BLOCKSIZE:(i+1)*BLOCKSIZE)=cb(:,:,counter);
            decodedB(j*BLOCKSIZE+1:(j+1)*BLOCKSIZE,1+i*BLOCKSIZE:(i+1)*BLOCKSIZE)=cr(:,:,counter);
            counter=counter+1;
        end
    end
%% SKLADANIE OBRAZU
img(:,:,1)=(decodedR);
img(:,:,2)=(decodedG);
img(:,:,3)=(decodedB);
%% WYSWIETLANIE 

img=uint8(img);
img= ycbcr2rgb(img);
I = ycbcr2rgb(I);
err = immse(I,img);
imwrite(img,'Po_kompresji.jpg');
montage({I,img});

%% FUNKCJE

function C = zigzag(A)
r = size(A,1); % number of rows of the input matrix
c = size(A,2); % number of columns of the input matrix
kk = 2;
C = [];
while kk <= r+c % the lowermost diagonal has only one element
                % which has the index r by c
                % scan the matrix as long as the last
                % diagonal is retrieved
    B = [];
    % iterate through every element of the input matrix
for ii = 1:r
    for jj = 1:c
        if ii + jj == kk % sum of the indices in a particular diagonal
                         % are the same
            B = [B,A(ii,jj)];
        end
    end
end
if mod(kk,2) == 0
    C = [C,flip(B)]; % reverse the order of the diagonal entries
                     % evenly
else
    C = [C,B];
end
kk = kk+1;
end
        
end



function [A] = izigsc(B,dim)
v = ones(1,dim); k = 1;
A = zeros(dim,dim);
for i = 1:2*dim-1
    C1 = diag(v,dim-i);
    C2 = flip(C1(1:dim,1:dim),2);
    C3 = B(k:k+sum(C2(:))-1);
    k = k + sum(C2(:));
    if mod(i,2) == 0
       C3 = flip(C3);
    end
    C4 = zeros(1,dim-size(C3,2));
    if i >= dim
       C5 = cat(2,C4, C3); 
    else       
        C5 = cat(2,C3,C4);
    end
    C6 = C2*diag(C5);
    A = C6 + A;
end
end

