extern	void	jostle_env(char*);

extern	void	jostle_wrkspc_input(int*,char*);

extern	void	jostle(int*,int*,int*,int*,int*,
		 int*,int*,int*,int*,int*,int*,int*,double*);

extern	void	pjostle_init(int*,int*);
extern	void	pjostle(int*,int*,int*,int*,int*,int*,int*,int*,
		 int*,int*,int*,int*,int*,int*,int*,double*);

extern	int	jostle_mem(void);
extern	int	jostle_cut(void);
extern	double	jostle_bal(void);
extern	double	jostle_tim(void);

