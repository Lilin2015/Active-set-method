function t = ASM( Params )
    %% prepare
    % read from Params
    W = Params.W;
    L = Params.L;
    b = Params.b;
    tw= Params.tw;
    
    % small number
    zero_t = 0.00000001;
    
    % pre-calc
    G = W + L;
    p_num = size(b,1);
    
    % initial working set and starting point can be changed, here are two examples
    % note that initial working set must be the active constraint at the starting point
    % example A
    t = G\(W*tw);         % the optimal solution without constraints
    w_set = find(t<=b);   % initial work set
    t = max(b,t);         % starting point
    % example B
%     t = tw;               % this works only if tw>=b, it is true in weighted dark channel dehazing
%     w_set = find(t<=b);
    
    %% loop
    iter = 0;
    while 1
        iter = iter + 1;
        w_num = size(w_set,1); 
        A = sparse(1:1:w_num,w_set,ones(w_num,1),w_num,size(tw,1));
        g = 2*(W*(t-tw)+L*t);
        B = [2*G,-A';A,zeros(w_num,w_num)];
        C = [-g;zeros(w_num,1)];
%         the function
%         B*PL=C
%         [2G -A'][     p]=[-g]
%         [ A  0 ][lambda] [ 0]
        PL= B\C;    p = PL(1:p_num);    lambda = PL(p_num+1:end);  % solve p and lambda 
        p(abs(p)<zero_t)=0;             % small numbers turn into zeros
        lambda(abs(lambda)<zero_t)=0;
        
        if all(~p) 
            % p is a null vector
            if ~all(lambda>=0)
                % all entries of lambda is non-negative
                minL = find(lambda==min(lambda));
                w_set(minL(1))=[];  % remove a constraint from the working set
                fprintf('constraint -,\t iter=%d\t workset=%d\t energy=%0.2f\n',iter, size(w_set,1), (t-tw)'*W*(t-tw)+t'*L*t);
                continue;
            else
                fprintf('problem solved,\t iter=%d\t workset=%d\t energy=%0.2f\n',iter, size(w_set,1), (t-tw)'*W*(t-tw)+t'*L*t);
                break  % all entries of lambda are non-negative, problem solved                
            end
        else
            % p is not a null vector
            blockRef = (b-t)./p;
            blockRef(w_set) = 2;
            blockRef(p>=0) = 2;
            minBR = min(blockRef);
            minBRF = find(blockRef == minBR);   % the indexes of blocking contraints
            alpha = min(1,minBR);   % calc alpha, the step gain
            t = t + alpha*p;
            if minBR >= 1
                % no blocking constraint, continue
                    fprintf('constraint h,\t iter=%d\t workset=%d\t energy=%0.2f\n',iter, size(w_set,1), (t-tw)'*W*(t-tw)+t'*L*t);
                continue;
            else
                % add a blocking constraint into working set
                w_set = [w_set;minBRF(1)];
                    fprintf('constraint +,\t iter=%d\t workset=%d\t energy=%0.2f\n',iter, size(w_set,1), (t-tw)'*W*(t-tw)+t'*L*t);
                continue;
            end
        end
    end
end

