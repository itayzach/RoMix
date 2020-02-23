
function L = laplacian(DATA, TYPE, options)  

% Calculate the graph laplacian of the adjacency graph of data set DATA.
%
% L = laplacian(DATA, TYPE, PARAM)  
% 
% DATA - NxK matrix. Data points are rows. 
% TYPE - string 'nn' or string 'epsballs'
% options - Data structure containing the following fields
% NN - integer if TYPE='nn' (number of nearest neighbors), 
%       or size of 'epsballs'
% 
% DISTANCEFUNCTION - distance function used to make the graph
% WEIGHTTYPPE='binary' | 'distance' | 'heat'
% WEIGHTPARAM= width for heat kernel
% NORMALIZE= 0 | 1 whether to return normalized graph laplacian or not 
%
% Returns: L, sparse symmetric NxN matrix 
%
% Author: 
%
% Mikhail Belkin 
% misha@math.uchicago.edu
%
% Modified by: Vikas Sindhwani (vikass@cs.uchicago.edu)
% June 2004

disp('Computing Graph Laplacian.');


NN=options.NN, 
DISTANCEFUNCTION=options.GraphDistanceFunction;
WEIGHTTYPE=options.GraphWeights;
WEIGHTPARAM=options.GraphWeightParam;
NORMALIZE=options.GraphNormalize;




% calculate the adjacency matrix for DATA
if strcmp(TYPE, 'nn')
    A = adjacency(DATA, TYPE, NN, DISTANCEFUNCTION);
elseif strcmp(TYPE, 'kernel')
    A = calckernel(options.Kernel,options.KernelParam,DATA);
else
    error('invalid option')
end
W = A;

% disassemble the sparse matrix
[A_i, A_j, A_v] = find(A);

switch WEIGHTTYPE
    
case 'distance'
   for i = 1: size(A_i)  
       W(A_i(i), A_j(i)) = A_v(i);
   end;
  
case 'binary'
 disp('Laplacian : Using Binary weights ');
    for i = 1: size(A_i)  
       W(A_i(i), A_j(i)) = 1;
    end;
 
case 'heat' 
    disp(['Laplacian : Using Heat Kernel sigma : ' num2str(WEIGHTPARAM)]);
    t=WEIGHTPARAM;
    for i = 1: size(A_i)  
       W(A_i(i), A_j(i)) = exp(-A_v(i)^2/(2*t*t));
    end;
case 'my_heat'
    fprintf('*******\n')
    fprintf('My heat means W = A = K which is the gaussian kernel gram matrix\n');    
    fprintf('*******\n')
    
%     [mPhi_W, Lambda_W] = eig(W);
                    
                    
%     fig1 = figure(2);
%     x1 = DATA(:,1);
%     x2 = DATA(:,2);
%     [mX1, mX2] = meshgrid(x1, x2);
%     for m = 0:3        
%         [vPhi_m_x1, ~] = SqExpEig(options.a_k, options.b_k, m, DATA(:,1));
%         [vPhi_m_x2, ~] = SqExpEig(options.a_k, options.b_k, m, DATA(:,2));
%         
%         % outter product since phi(x1,x2)=phi(x1)phi(x2)
%         vPhi_m_x1x2 = vPhi_m_x1 .* vPhi_m_x2; 
%         
%         vPhi_W = mPhi_W(:,m+1);
%         
%         subplot(2,2,m+1);
% %         surf(mX1, mX2, mPhi_m_x1x2, 'edgecolor', 'none')
% %         hold on
%         scatter3(DATA(:,1), DATA(:,2), vPhi_m_x1x2, 'filled', 'bo')
%         hold on
%         scatter3(DATA(:,1), DATA(:,2), vPhi_W, 'filled', 'ro')
% %         view(2)
% %         colorbar()
% %         xlabel('$x_1$', 'Interpreter', 'latex')
% %         ylabel('$x_2$', 'Interpreter', 'latex')
% %         zlabel(['$\phi_' num2str(m) '(x_1,x_2)$'], 'Interpreter', 'latex')
%     end 
        
otherwise
    error('Unknown Weighttype');   
end

D = sum(W(:,:),2);   

if NORMALIZE==0
    L = spdiags(D,0,speye(size(W,1)))-W;
else % normalized laplacian
    D=diag(sqrt(1./D));
    L=eye(size(W,1))-D*W*D;
end
