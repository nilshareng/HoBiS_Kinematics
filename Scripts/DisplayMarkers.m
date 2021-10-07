function DisplayMarkers(markers,N)

figure(N);
hold on;

if size(fieldnames(markers),1)==12
    
    plot3(markers.RFWT(1,1), markers.RFWT(1,2), markers.RFWT(1,3),'bx');
    plot3(markers.LFWT(1,1), markers.LFWT(1,2), markers.LFWT(1,3),'bx');
    plot3(markers.RBWT(1,1), markers.RBWT(1,2), markers.RBWT(1,3),'bx');
    plot3(markers.LBWT(1,1), markers.LBWT(1,2), markers.LBWT(1,3),'bx');
    plot3(markers.LKNE(1,1), markers.LKNE(1,2), markers.LKNE(1,3),'bx');
    plot3(markers.RKNE(1,1), markers.RKNE(1,2), markers.RKNE(1,3),'bx');
    plot3(markers.LKNI(1,1), markers.LKNI(1,2), markers.LKNI(1,3),'bx');
    plot3(markers.RKNI(1,1), markers.RKNI(1,2), markers.RKNI(1,3),'bx');
    plot3(markers.LANE(1,1), markers.LANE(1,2), markers.LANE(1,3),'bx');
    plot3(markers.RANE(1,1), markers.RANE(1,2), markers.RANE(1,3),'bx');
    plot3(markers.LANI(1,1), markers.LANI(1,2), markers.LANI(1,3),'bx');
    plot3(markers.RANI(1,1), markers.RANI(1,2), markers.RANI(1,3),'bx');
    
    xlim([-300 300]);
    ylim([-300 300]);
    zlim([-1000 200]);
    
elseif size(fieldnames(markers),1)==24
    
    plot3(markers.RFWT(1,1), markers.RFWT(1,2), markers.RFWT(1,3),'bx');
    plot3(markers.LFWT(1,1), markers.LFWT(1,2), markers.LFWT(1,3),'bx');
    plot3(markers.RBWT(1,1), markers.RBWT(1,2), markers.RBWT(1,3),'bx');
    plot3(markers.LBWT(1,1), markers.LBWT(1,2), markers.LBWT(1,3),'bx');
    plot3(markers.LKNE(1,1), markers.LKNE(1,2), markers.LKNE(1,3),'bx');
    plot3(markers.RKNE(1,1), markers.RKNE(1,2), markers.RKNE(1,3),'bx');
    plot3(markers.LKNI(1,1), markers.LKNI(1,2), markers.LKNI(1,3),'bx');
    plot3(markers.RKNI(1,1), markers.RKNI(1,2), markers.RKNI(1,3),'bx');
    plot3(markers.LANE(1,1), markers.LANE(1,2), markers.LANE(1,3),'bx');
    plot3(markers.RANE(1,1), markers.RANE(1,2), markers.RANE(1,3),'bx');
    plot3(markers.LANI(1,1), markers.LANI(1,2), markers.LANI(1,3),'bx');
    plot3(markers.RANI(1,1), markers.RANI(1,2), markers.RANI(1,3),'bx');
    plot3(markers.RFem9(1,1),markers.RFem9(1,2),markers.RFem9(1,3),'bx');
    plot3(markers.LFem9(1,1),markers.LFem9(1,2),markers.LFem9(1,3),'bx');
    plot3(markers.RFem6(1,1),markers.RFem6(1,2),markers.RFem6(1,3),'bx');
    plot3(markers.LFem6(1,1),markers.LFem6(1,2),markers.LFem6(1,3),'bx');
    plot3(markers.RTib6(1,1),markers.RTib6(1,2),markers.RTib6(1,3),'bx');
    plot3(markers.LTib6(1,1),markers.LTib6(1,2),markers.LTib6(1,3),'bx');
    plot3(markers.RTib1(1,1),markers.RTib1(1,2),markers.RTib1(1,3),'bx');
    plot3(markers.LTib1(1,1),markers.LTib1(1,2),markers.LTib1(1,3),'bx');
    plot3(markers.RTal1(1,1),markers.RTal1(1,2),markers.RTal1(1,3),'bx');
    plot3(markers.LTal1(1,1),markers.LTal1(1,2),markers.LTal1(1,3),'bx');
    
    xlim([-300 300]);
    ylim([-300 300]);
    zlim([-1000 200]);
    
end

% % end
% [markers.RFWT(1,1), markers.RFWT(1,2), markers.RFWT(1,3);
% markers.LFWT(1,1), markers.LFWT(1,2), markers.LFWT(1,3);
% markers.RBWT(1,1), markers.RBWT(1,2), markers.RBWT(1,3);
% markers.LBWT(1,1), markers.LBWT(1,2), markers.LBWT(1,3);
% markers.LKNE(1,1), markers.LKNE(1,2), markers.LKNE(1,3);
% markers.RKNE(1,1), markers.RKNE(1,2), markers.RKNE(1,3);
% markers.LKNI(1,1), markers.LKNI(1,2), markers.LKNI(1,3);
% markers.RKNI(1,1), markers.RKNI(1,2), markers.RKNI(1,3);
% markers.LANE(1,1), markers.LANE(1,2), markers.LANE(1,3);
% markers.RANE(1,1), markers.RANE(1,2), markers.RANE(1,3);
% markers.LANI(1,1), markers.LANI(1,2), markers.LANI(1,3);
% markers.RANI(1,1), markers.RANI(1,2), markers.RANI(1,3);
% markers.RFem9(1,1),markers.RFem9(1,2),markers.RFem9(1,3);
% markers.LFem9(1,1),markers.LFem9(1,2),markers.LFem9(1,3);
% markers.RFem6(1,1),markers.RFem6(1,2),markers.RFem6(1,3);
% markers.LFem6(1,1),markers.LFem6(1,2),markers.LFem6(1,3);
% markers.RTib6(1,1),markers.RTib6(1,2),markers.RTib6(1,3);
% markers.LTib6(1,1),markers.LTib6(1,2),markers.LTib6(1,3);
% markers.RTib1(1,1),markers.RTib1(1,2),markers.RTib1(1,3);
% markers.LTib1(1,1),markers.LTib1(1,2),markers.LTib1(1,3);
% markers.RTal1(1,1),markers.RTal1(1,2),markers.RTal1(1,3);
% markers.LTal1(1,1),markers.LTal1(1,2),markers.LTal1(1,3)];



