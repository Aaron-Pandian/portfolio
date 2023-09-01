import pygame
pygame.init()

win = pygame.display.set_mode((500,500))
pygame.display.set_caption('My First Game')

x = 40
y = 40
w = 40
h = 40
vel = 5

run = True
while run:
    pygame.time.delay(25)
    win.fill((0, 0, 0))
    pygame.draw.rect(win, (250, 0, 0), (x,y,w, h))
    pygame.display.update()

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            run = False

    keys = pygame.key.get_pressed()
    if keys[pygame.K_LEFT] and x > 0:
        x -= vel
    if keys[pygame.K_RIGHT] and x < 500 - w:
        x += vel
    if keys[pygame.K_UP] and y > 0:
        y -= vel
    if keys[pygame.K_DOWN] and y < 500 - h:
        y += vel

pygame.quit()
