CDF      
      
len_string     !   len_line   Q   four      	time_step          
num_qa_rec        num_info   o   len_name   !   num_dim       	num_nodes         num_elem      
num_el_blk        num_side_sets         num_el_in_blk1        num_nod_per_el1       num_side_ss1      num_side_ss2      num_side_ss3      num_side_ss4      num_nod_var          	   api_version       @�ff   version       @�ff   floating_point_word_size            	file_size               int64_status             title         Default Sierra Title       processor_info                   last_written_time         ?�ߐU�|   maximum_name_length          
         
time_whole                            1�   
qa_records                             �      �   info_records                      #       	   node_num_map                    d      ,$   elem_num_map      	              @      ,�   	eb_status         
                    ,�   eb_prop1      
         name      ID              ,�   eb_names      
                 $      ,�   	ss_status                             ,�   ss_prop1               name      ID              -   ss_names                       �      -   coordx                      �      -�   coordy                      �      .`   
coor_names                         D      /(   connect1                  	elem_type         QUAD4               /l   elem_ss1                          0l   side_ss1                          0|   elem_ss2                          0�   side_ss2                          0�   elem_ss3                          0�   side_ss3                          0�   elem_ss4                          0�   side_ss4                          0�   vals_nod_var1                          �      1�   vals_nod_var2                          �      2d   vals_nod_var3                          �      3,   vals_nod_var4                          �      3�   vals_nod_var5                          �      4�   name_nod_var                       �      0�Conchas                          4.25.4-137-gfb9a5240             2012/05/16                       05:53:46                         Node: wsblade007, OS: Linux 2.6.18-238.9.1.el5, #1 SMP Fri Mar 18 12:42:39 EDT 2 Ioex_DatabaseIO.C 2012/04/19 gdsjaar                                             /scratch/bcarnes/SIERRA_TRILINOS/code/results/2012_05_16_05_52_16/conchas_rtest_ $ Aprepro (Revision: 2.26) Wed May 16 05:53:46 2012                              Begin Sierra Conchas                                                               Begin Finite Element Model femodel                                                 Database Name = input1.exo                                                       Database Type = ExodusII                                                                                                                                          Begin Parameters For Block block_1                                                 Material material_air                                                          End                                                                            End                                                                                                                                                               Begin Property Specification For Conchas Material material_air                     gamma = 1.4                                                                      specific_r = 287.097384766765                                                  End                                                                                                                                                               load user plugin file ./mms_user.so                                                                                                                               Begin Conchas Procedure conchas_procedure                                           Begin Solution Control Description                                                  Use System Main                                                                  Begin System Main                                                                   Begin Transient The_Transient_Block                                                 Advance conchas_region                                                        End                                                                           End system main                                                                                                                                                   Begin Parameters For Transient The_Transient_Block                                  Start Time = 0.0                                                                 Termination Time = 100.0                                                         number of steps = 5000                                                           BEGIN PARAMETERS FOR CONCHAS REGION conchas_region                                  TRANSIENT STEP TYPE IS automatic                                                 cfl = 2.0                                                                        cfl ramp increment = 0.01 max = 100                                           END                                                                           END                                                                           end                                                                                                                                                              Begin Conchas Region conchas_region                                                                                                                                 Use Finite Element Model femodel                                                                                                                                  Activate Edge Based Algorithm                                                                                                                                     Begin Solution Options                                                              Activate Equation Euler                                                                                                                                           Number of PI Sweeps = 5                                                                                                                                           Upwind Method = MUSCL For Equation Euler                                                                                                                          Upwind Limiter = Barth For Equation Euler                                                                                                                         Freeze Limiter at Step 500                                                                                                                                        Minimum Number Of Nonlinear Iterations = 1                                       Maximum Number Of Nonlinear Iterations = 1                                                                                                                        COMPUTE STEADY SOLUTION USING PSEUDO TRANSIENT METHOD                            NONLINEAR RESIDUAL NORM TOLERANCE = 5.0e-8                                                                                                                        Activate Manufactured Solution Source Terms                                      #Load Exact Pressure euler_p From File mms_euler.so                           End                                                                                                                                                               Begin Flow State freeStream                                                        pressure = 35651.28116                                                           temperature = 236.215                                                            direction = 1, 0, 0                                                              mach number = 2.5                                                              End                                                                                                                                                               Begin Initial Condition Block icblock                                              Volume is block_1                                                                Use Flow State freeStream                                                      End                                                                                                                                                               # inflow (left)                                                                  Begin Supersonic Inflow on surface surface_1                                       Use Flow State freeStream                                                      End                                                                                                                                                               # top wall                                                                       Begin Symmetry Boundary Condition on surface surface_2                           End                                                                                                                                                               # bottom wall                                                                    Begin Symmetry Boundary Condition on surface surface_4                           End                                                                                                                                                               # outflow (right)                                                                Begin Open Boundary Condition on surface surface_3                               End                                                                                                                                                               Begin Results Output Label output                                                  Database Name = output1.e                                                        At Step 0, Increment = 100                                                                                                                                        Nodal variables = density                                                        Nodal variables = velocity                                                       Nodal Variables = pressure as P                                                  Nodal Variables = temperature as T                                             End                                                                                                                                                             End                                                                            End                                                                            End                                                                                                          	   
                                                                        	   
                        block_1                                                     surface_1                        surface_2                        surface_3                        surface_4                        @       @       ?�      ?�      @       ?�      @       ?�      @       ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�                                                      ?�      ?�              ?�      ?�      ?�      ?�      ?�      ?�      ?�              ?�      ?�      ?�      ?�              ?�      ?�      ?�      ?�              ?�      ?�      ?�      x                                y                                                                            	   
                                             
                                                                                                                                    	                                                                  density                          P                                T                                velocity_x                       velocity_y                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          ?�<�H�o�?�����Q?�ذ��?�ՂI�B�?�Ө����?����߫A?᪰����?��b��`O?����Y�?�Т6�?�Z�x�?������?�	�1C}J?�E��x�/?�{Ơ�?��z�7T?ߏ3nLq�?��.���?�	��?���0}��?���8D�^?�/\H,�*?��_#��?�9Wj'�B?�U��k�Y?�.���٠@�R��_�@�N�OB@��o�&@�~�|��`@�o&��d�@�j�üY@�y}%�I@�X$�7�@�~�3�X�@�/���@⍯�%��@��"��&@�/����@�X�#eK@߉	J�@���tm�@�FZCի�@�:�'w�@����b@␁0�gs@�A�@ߚM����@�f^pC��@�|;J@�p*�6@ms��;>@m�M�k@l��@m�� m^@m~A���>@law=5]@m{�̉(�@lo�t�o@m�~�Z\�@mE�6�̍@ma��a|D@m�41���@m�oᖇ@lL<��R@m�G��1@n��6D@m>�"衔@n�MY�?@nl��1��@m����@m���|�(@l� Ur�O@n���@n�x��A�@m��%�G�@�ZТ9�@����62@���2�0�@���@��@��JV��@��j���@��K��@�T�W��@��g�f/@�-l~Ţ@�@�>r��@�X�6�N(@���-,�@��Gřn@�����V�@�8p�xGx@��!0�@��	"��@��eEV/T@��h L�@����F"^@�J���]@���W�@�w_����@�����&�?�$��b�h@l�y��#�3&4\ ?���$�Q:@S�v����9X�Q��@>O���m�2��K�*?��Xo�E��1.�����!�C��������]t��%�%�z�v�)l�K��(����I@59��N�3����!%W{@8�� SX,@1/���-?�郊��@0��\�����WD�@7��F�A@/��<������?��-���?�Ȓ�i�@?��ճ?���o�U?����'?��)Rd�?��D�~(?��i\��?����ܶV?�Ѝ�ή)?�Z��gk�?�8�.-.?��}��]?�E�:� �?�{�LS�??��qzz�?ߕ�|�a�?߱����?���i�?����>?����?�]?�.�Q��?� F�@<,?�9����?�U�\'�?�.�G��#@�Re��T@�OM�M�@��;Y���@�~�frM@�o`�[�\@ᒍI���@�y!��@�u�w>�@�~�P#�@�n7<:@⎡��[�@��f�4�@�	���z@�X�� �@@߉Ovc<@����@�:(@�;��E��@���3�@␨��@�@���U;@ߚ��XT�@�f��	��@�r�@�{�
�@mr���';@m�#���q@l�˘�@m�g����@m~R�&NO@la�=�g@m{����+@lo�e��@m�w@�g@mEs`�t`@mb?ʕ$p@m�/dN/�@mЋ;&@lL7�ӻ�@m��%g@n�Y��@m7�F��>@n����5I@nl�N?�@m��VN�@m�u-F�]@l�'4�
)@n�6}��+@n���x^@m���g@�E<��>@��^s�@�����8�@���{�\i@�ݧs��@��˿�@���t@�T��k��@��m5p@�.��·@�@ؓۅ@�X���	?@��v��e@����@����O�q@�7�(̣�@�C����@���|�U@��c�5hU@�떟qj�@�����F@�J2q�%�@���-|@�w\p�+%@���v�V�?��2�w@{��|!��3&A ?�I�9�Z@l*�G��9WYs=�@C��m��2�r��2?����gEq��3� Mg�!U�f������%�Lnw��)l��H���#�s�@5f����l�59�@8�VC<ʔ@1/�*�|r?�}Zʪg@0��$�f����5�ס@7�٩x{�@/w�t����[h�?�Kϟ�ź?���;'?���?���o��{?���z9v?�ݻ_r�?᪛��Ƶ?�鋻�Y�?�͹��l?�В�a�?�[���;?��M���?�-���+?�E/'8�?�{�D
2�?����?߃�9�Ю?���C�k2?���\?��S�<?�Л���?�/Q�3L?�#�9nb�?�9[��]?�U���0?�.�ɧ��@�R�n\G@�O#<K@���g��Y@�}�?���@�n��5}w@�O�ak@�y����@��ݕ�@�~���Ͷ@�Uo>$O@��s��.@��(��F!@�}|�G @�X�?�@߇���W�@��"��= @�Z�W�B�@�:#2܂�@�:���@�&s@�@�VYs�@ߠHk�5V@�f^}y>&@�Vj��@�iDa@ms#Y��@m�	"��@l��d��@m��o�(@m~'$B@lam����@m{�����@lo��)C@m���F�9@mE��k;6@ma��K,@m�]&��@m�j;�]@lK��0@mm2&'l@n} �v�@mH4Bx:6@n�?�}GI@nlmB<�@m�vrG��@m�	t��@l��;\��@n���&@n��{��@m��(p�@�[�{lU@��QѪc@�����@��޲1@��z@����L@��	��@�Tȡ~@���=<N@�6DɏU@�@��{ߜ@�X�I�l@����@��	���?@��£���@�8t\��@����0@�Ԥ�t�@��u?ܛ�@��C���@���4ye�@�J��kA8@�}P@4�@�wZ�}��@���3�?��l��@{��P��3&� K	V?��A��@F�/i@9�9V���@T���
��2�\e�U?��7v�~�� ���c`�!V�Yqx���=9�t�%�h�X��)x�;F�A�Hx���@4��㶇��R��<�@8��Wn��@1(P����?�J�ez@0�D&�Iؿ��3b��@7�M}�e�@/����5���;���?�p(!�Xa?����B�?൰���?�գ}�?����?��+	���?���t��?��C���?���,\��?�����;W?�YA�|�	?�<���9?��u?�F9B[C?�{��|�?��	�;W?ߐk
�?߾����?����?���r�?����\tI?�/�g3$?�ID]?�9L���I?�U�o���?�.�φ�@�R�#}�"@�N���{@��׵q��@�����@�oae/�/@ᒽ�9"�@�yE�9�@�ke��@�D~�G@�,Eu�g@�g:��<@��d$���@�tn��<@�X�l���@߈��@��kW���@�D�I%Y�@�:uWm�@��w�1u@␫�p�)@�A�`��@ߘ_eF@�fs a��@�,<�<�@����@mr�a͑%@m�L�sK�@l�2�KY@m����\@m~Px�s�@la��\J�@m{��*v�@lo�9@v�@m����fL@mE*��sK@ma܂�"	@m��cZ@m�y'!@lLP�w#y@m}��iC@n�jd;@m>�k�x�@n�ZVa�5@nlQ�u]5@m���vR�@m��sm�@l����=�@n�RƉt�@n�~'�6)@m��g�6:@�iЅΑ@��;��@@��x,��u@��x��@�ᘢD�@���F@���d1�@�T���M@���k�W@�%k�e@�@���f-@�X�#�>�@��9>&F@���iY�@����^��@�8�����@�&��M7@�}�W�P@��c���0@�뛈��^@���Zt\@�I�x��@�t�m�@�wY ��@����P�?���]���@�{�]aG�3!z0�ӊ?�R�b��	@�w��r��9S�˱@A��{�a�2ĕ��?�I޻y/������"�!��P����c��nn��%����$�)ua�����.J9@5Sr0�9u��EeX9�@8�:�"�@10�3���?���ۨ@0�ZLV���>�N�@7���@/�&������Xs԰�?Ķ��?��D���?�� Q�?���+#N?��g�2�?�ާLu�?� 44�P?��g!�]?��.��u�?��BՓ�?�Z��A��?����Cx?����g?�E��%>?�{�PBH?���G}?ߗY�S�?߭S���?�9:��?��2Ǣ�?��Xɬ?�-��Hn'?� �{Y�?�9�[��H?�U��F�?�.���k@�R�XL�@�O
��A�@���`8R@�X��7�@�p���5@���]�@�yx��A8@࿾��{U@�~��v�|@��I%�@�
1�Y@��j��*P@��!��@�XʦE��@߉-ÛI@�㨜�@@�6��|�@�<k*�@���:�@����|�@�@��	U@ߛi�b�@�f�6�ic@�2ҙ�O@���t
�@mr�E� �@m�(jva@l�&&: @m�=�'��@m~�Y^D�@la���^@m{����@lp(~O�^@m�zs}c�@mEHrJ-@mbf\+5d@m����@m�X���@lLX�!�@mКN�o@n�r�b$@m63��DF@n��`�y@nl��B��@m��D�t�@m�Y4{�Z@l���Zk�@n�`6N��@n�����(@m��ݛ��@�:}n��@��<��@�����n@���e�<@��#�8@�묖�^@��+�-@�T�1W�k@�����y@�9�l�@�@�s��@�X\ZP�@����"�@��.��@������@�7��(@@��K�^�@���H+�@��b�@�[@��X�V�@���|�I@�JQ:��V@�{�h�$@�wW���@������?�}�e;�@��em���3$���H�?�Im�?�g@����X�9TӾQ@W��	*�2��wM?�{?Z"o��72��u�!H�]ο��D܈�%����Ys�)qJ�PBV����Ǡ@5d�+�������@8�k�"��@10VN���?���G,�@0�b��y����A�:�@7��̀L@/�5������;�b?�ߐU�|?���N�e?��P��?��&��+A?���0>?����a�?�<� ��?��8��;N?��VT��_?���C	+?�H/��Y�?�z���?��F9F�?�C~�+?�u����?���te�x?߉�]D?� \���?���b��?��GeҦ?��J����?��	l��?�"�y< ?�:`����?�Y���?�3�Go��@�B�I@�X�@���U�@�i���.@�r�lM�"@�1��&C@�w0�:@��k�a@�x�]xU6@�2@]n@�8U��@�БGn�*@ፓ�^�}@�ONH���@�p���s�@�۞��@�ok�Z��@�2٬�)@����{@��M/�@�2�c��@ߦ"�H@�gyv�C@�āk�h@�-��f�@mk�	��@m����c@l��G��b@m���W�@m~<�&@la���Z�@my���b@ll��א�@m�;O��@m?I���I@mkx�|	@m��Bz��@m$�@lFSb�x@m�A<o@n�MAW�b@mO,G�0�@n�.T�z�@njM�{�K@m��|c�@m��_�{p@m ��,��@n���:�@n�ߞ�9@m𦋎Z@���6�@��\j@����n@��'��eo@�B�.g@�qh?�@���F��@�U�>b�r@��ݡ�Wo@���F��@�=�t֥�@�X,���@�&�n @��{ޝV@��FY1�@�8��ѵ#@�h�<@��9�k@��S��@��\۽@��"�`#6@�F�Yu��@���0K�@�u�\��@��^�1?�y�#� @�"�l�x�2���=&@jDb�C}@�#Q �8��J|@f]l���2�~-|\�?��ϳ�����Ϝ�
e� 䲊x(@��}��%��^��*'���ʪ��%l7��@4ͺ��UE@>&gږ�@8B�V��@0��#0?��Ё��i@0����V�,��@7��E�R�@/h��"ܿ���͸