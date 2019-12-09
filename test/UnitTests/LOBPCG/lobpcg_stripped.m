function [X,lambda,varargout] = lobpcg_stripped(X,A, tol, varargin)
    [n,blockSize]=size(X);
    verbosityLevel = 1;
    failureFlag = 1;
    
    if nargin < 4
        maxIterations = min(n, 20);
    else
        maxIterations = varargin{1};
    end

    %Make X orthonormal
    % [X, flag] = cholQ(X);
    % assert(flag == 0, 'The initial approximation is not full rank')
    
    % Preallocation
    condestGhistory=zeros(1,maxIterations+1);
    P=zeros(n,blockSize, 'like', X);
    AP=zeros(n,blockSize, 'like', X);
    
    [X, lambda, AX, R] = RR(X, A);
    assert(norm(X'*X - eye(size(X, 2)), 'fro') <= 1e-8)
    condestGhistory(1)=-log10(eps)/2;  %if too small cause unnecessary restarts
    active_indices = true(blockSize,1);
    
    for iterationNumber=1:maxIterations
        residualNorms = vecnorm(R)';
        
        active_indices = (residualNorms > tol) & active_indices;
        % sum(active_indices)
        R = R(:, active_indices);
        P = P(:, active_indices);
        AP = AP(:, active_indices);
        
        if sum(active_indices) == 0
            failureFlag=0; %all eigenpairs converged
            break
        end
        % Making active residuals orthogonal to X
        % assert(norm(X*(X'*R)) < 1e-8)
        R = R - X*(X'*R);
        
        %Making active residuals orthonormal
        [R, flag] = cholQ(R);
        if flag ~= 0
            warning('BLOPEX:lobpcg:ResidualNotFullRank',...
                'The residual is not full rank.');
            break
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [P, AP, X, AX, lambda, condestGhistory] = RayeleighRitz(A, ...
            X, AX, lambda, R, P, AP, residualNorms, condestGhistory, ...
            iterationNumber, blockSize);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%end RR
        
        if false % verbosityLevel
            fprintf('Iteration %i current block size %i \n',...
                iterationNumber, sum(active_indices));
            fprintf('Eigenvalues lambda %17.16e \n',lambda);
            fprintf('Residual Norms %e \n',residualNorms');
        end

        if iterationNumber == maxIterations
            fprintf('Warning: maximum iterations reached without converging\n')
        else
            R = AX - X.*lambda';
        end
    end
    %The main step of the method was the CG cycle: end

    %Postprocessing
    %Making sure X's are "exactly" othonormalized by final "exact" RR
    [X, lambda, ~, R] = RR(X, A);
    residualNorms=vecnorm(R)';

    varargout(1)={failureFlag};
    varargout(2)={residualNorms};
end


function [P, AP, X, AX, lambda, condestGhistory] = RayeleighRitz(A,...
    X, AX, lambda, R, P, AP, residualNorms, condestGhistory,...
    iterationNumber, blockSize)

    AR = A*R;
    
    % The Raileight-Ritz method for [X R P]
    if residualNorms <= eps^0.6  %suggested by Garrett Moran, private
        explicitGramFlag = 1;
    else
        explicitGramFlag = 0;
    end

    if iterationNumber == 1
        restart = 1;
    else
        %Making active conjugate directions orthonormal
        [P, AP, restart] = cholQ2(P, AP);
        if restart
            warning('BLOPEX:lobpcg:DirectionNotFullRank',...
                'The direction matrix is not full rank.');
        end
    end

    for cond_try=1:2           %cond_try == 2 when restart
        [oA, oB] = getGram(X, R, P, AX, AR, AP, lambda,...
            blockSize, explicitGramFlag, restart);
        
        condestG = log10(cond(oB))+1;
        condestGmean = mean(condestGhistory(max(1,iterationNumber-10-...
            round(log(size(R,2)))):iterationNumber));
        if (condestG/condestGmean > 2 && condestG > 2 )|| condestG > 8
            %black magic - need to guess the restart
            if true %verbosityLevel
                fprintf('Restart on step %i as condestG %5.4e \n',...
                    iterationNumber,condestG);
            end
            if cond_try == 1 && ~restart
                restart=1; %steepest descent restart for stability
            else
                warning('BLOPEX:lobpcg:IllConditioning',...
                    'Gramm matrix ill-conditioned: results unpredictable');
            end
        else
            break
        end
    end
    [W, L]=eig(oA,oB);
    % Shouldn't we sort the eigenvalues??
    % assert(norm(sort(diag(L)) - diag(L)) == 0)
    W = W(:, 1:blockSize);
    lambda=diag(L(1:blockSize,1:blockSize));
    
    if ~restart
        % X = [X R P]*W;
        P = [R P]*W(blockSize+1:end, :);
        AP = [AR AP]*W(blockSize+1:end, :);
    else %use block steepest descent
        % X = [X R]*W;
        P = R*W(blockSize+1:end, :);
        AP = AR*W(blockSize+1:end, :);
    end
   
    X = X*W(1:blockSize,:) + P;
    AX = AX*W(1:blockSize,:) + AP;
    clear W;
    %%end RR
    condestGhistory(iterationNumber+1) = condestG;
end

function [A, flag] = cholQ(A)
    [C, flag] = chol(A'*A);
    if flag == 0
        A = A/C;
    end
end

function [A, B, flag] = cholQ2(A, B)
    [C, flag] = chol(A'*A);
    if flag == 0
        A = A/C;
        B = B/C;
    end
end

function [X, l, AX, R] = RR(X, A)
    XX=X'*X;
    XX=(XX' + XX)*0.5;
    AX = A*X;
    XAX = X'*AX;
    XAX = (XAX + XAX')*0.5;
    [V, L] = eig(XAX,XX);
    l=diag(L);
    X = X*V;
    AX = AX*V;
    R = AX - X.*l';
end

function [oA, oB] = getGram(X, R, P, AX, AR, AP, lambda,...
    blockSize, explicitGramFlag, restart)

    oXAR=AX'*R;
    oRAR=AR'*R;
    oRAR=(oRAR'+oRAR)*0.5;
    
    if ~restart
        oXAP=AX'*P;
        oRAP=AR'*P;
        oPAP=AP'*P;
        oPAP=(oPAP'+oPAP)*0.5;
        oXP=X'*P;
        oRP=R'*P;
    end

    if explicitGramFlag
        oXAX=AX'*X;
        oXAX=(oXAX'+oXAX)*0.5;
        oXX=X'*X;
        oXX=(oXX'+oXX)*0.5;
        oRR=R'*R;
        oXR=X'*R;
        oRR=(oRR'+oRR)*0.5;

        if ~restart
            oA = [oXAX      oXAR     oXAP
                  oXAR'     oRAR     oRAP
                  oXAP'     oRAP'    oPAP];
        
            oPP=P'*P;
            oPP=(oPP'+oPP)*0.5;
            oB = [oXX    oXR    oXP
                  oXR'   oRR    oRP
                  oXP'   oRP'   oPP];
        else
            oA = [oXAX     oXAR
                  oXAR'    oRAR];
              
            oB = [oXX      oXR
                  oXR'     eye(size(R, 2))];
        end
    else
        if ~restart
            oA = [diag(lambda)   oXAR     oXAP
                  oXAR'          oRAR     oRAP
                  oXAP'          oRAP'    oPAP];
            oB = [eye(blockSize)                  zeros(blockSize,size(R, 2))   oXP
                  zeros(blockSize,size(R, 2))'    eye(size(R, 2))               oRP
                  oXP'                            oRP'                          eye(size(P, 2)) ];
        else
            oA = [diag(lambda)   oXAR
                  oXAR'          oRAR];
            oB = eye(blockSize+size(R, 2));
        end
    end
end