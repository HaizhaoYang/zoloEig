function [K, H, Q, Z] = util_hh2th(K, H)
% UTIL_HH2TH    Transform rational pencil to a polynomial one.
%
% The function brings a (quasi-)upper-Hessenberg pencil to an
% upper-triangular--upper-Hessenberg one.
%
% [K, H, Q, Z] = util_hh2th(K, H) finds orthogonal matrices Q and Z
% such that Q*K*Z is upper-triangular and Q*H*Z is
% upper-Hessenberg.  

  m = size(H, 2);
  Q = eye(m+1);
  Z = eye(m);

  j = 1;
  while j <= m
    if j+1 <= m && H(j+2, j) ~= 0

      % Annihilate element K(j+1, j) from the left.
      [s,c] = util_givens(K(j, j), K(j+1, j));
      K(j:j+1, :) = [c -s; s' c]*K(j:j+1, :);
      H(j:j+1, :) = [c -s; s' c]*H(j:j+1, :);
      Q(j:j+1, :) = [c -s; s' c]*Q(j:j+1, :);

      % Annihilate element K(j+2, j+1) from the left.
      [s,c] = util_givens(K(j+1, j+1), K(j+2, j+1));
      K(j+1:j+2, :) = [c -s; s' c]*K(j+1:j+2, :);
      H(j+1:j+2, :) = [c -s; s' c]*H(j+1:j+2, :);
      Q(j+1:j+2, :) = [c -s; s' c]*Q(j+1:j+2, :);

      for ii = j+2:-1:3
        if ii ~= 3
          % Annihilate element H(ii, ii-3) from the right.
          [s,c] = util_givens(H(ii, ii-2), H(ii, ii-3));
          H(:, ii-3:ii-2) = H(:, ii-3:ii-2)*[c -s; s' c];
          K(:, ii-3:ii-2) = K(:, ii-3:ii-2)*[c -s; s' c];
          Z(:, ii-3:ii-2) = Z(:, ii-3:ii-2)*[c -s; s' c];
        end
        % Annihilate element H(ii, ii-2) from the right.
        [s,c] = util_givens(H(ii, ii-1), H(ii, ii-2));
        H(:, ii-2:ii-1) = H(:, ii-2:ii-1)*[c -s; s' c];
        K(:, ii-2:ii-1) = K(:, ii-2:ii-1)*[c -s; s' c];
        Z(:, ii-2:ii-1) = Z(:, ii-2:ii-1)*[c -s; s' c];

        if ii ~= 3
          % Annihilate element K(ii-2, ii-3) from the left.
          [s,c] = util_givens(K(ii-3, ii-3), K(ii-2, ii-3));
          K(ii-3:ii-2, :) = [c -s; s' c]*K(ii-3:ii-2, :);
          H(ii-3:ii-2, :) = [c -s; s' c]*H(ii-3:ii-2, :);
          Q(ii-3:ii-2, :) = [c -s; s' c]*Q(ii-3:ii-2, :);
        end

        % Annihilate element K(ii-1, ii-2) from the left.
        [s,c] = util_givens(K(ii-2, ii-2), K(ii-1, ii-2));
        K(ii-2:ii-1, :) = [c -s; s' c]*K(ii-2:ii-1, :);
        H(ii-2:ii-1, :) = [c -s; s' c]*H(ii-2:ii-1, :);
        Q(ii-2:ii-1, :) = [c -s; s' c]*Q(ii-2:ii-1, :);
      end % ii = j+2:-1:3

    else

      % Annihilate element K(j+1, j) from the left.
      [s,c] = util_givens(K(j, j), K(j+1, j));
      K(j:j+1, :) = [c -s; s' c]*K(j:j+1, :);
      H(j:j+1, :) = [c -s; s' c]*H(j:j+1, :);
      Q(j:j+1, :) = [c -s; s' c]*Q(j:j+1, :);
      % Bulge-chasing part.
      for ii = j+1:-1:3
        % Annihilate element H(ii, ii-2) from the right.
        [s,c] = util_givens(H(ii, ii-1), H(ii, ii-2));
        H(:, ii-2:ii-1) = H(:, ii-2:ii-1)*[c -s; s' c];
        K(:, ii-2:ii-1) = K(:, ii-2:ii-1)*[c -s; s' c];
        Z(:, ii-2:ii-1) = Z(:, ii-2:ii-1)*[c -s; s' c];

        % Annihilate element K(ii-1, ii-2) from the left.
        [s,c] = util_givens(K(ii-2, ii-2), K(ii-1, ii-2));
        K(ii-2:ii-1, :) = [c -s; s' c]*K(ii-2:ii-1, :);
        H(ii-2:ii-1, :) = [c -s; s' c]*H(ii-2:ii-1, :);
        Q(ii-2:ii-1, :) = [c -s; s' c]*Q(ii-2:ii-1, :);
      end % ii = j+1:-1:3
    end
    j = j+1;
  end % while j <= m

end
